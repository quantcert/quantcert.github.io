/*
 Copyright (C) Tommi Junttila
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License version 2
 as published by the Free Software Foundation.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <cstdio>
#include <cstring>
#include <cstdarg>
#include <list>
#include "defs.hh"
#include "bc.hh"
#include "timer.hh"



#ifndef BC_HAS_MINISAT
/*
 *
 * No MiniSAT, a call will result in a terminating internal error.
 *
 */
int BC::minisat_solve(const bool perform_simplifications
		      , const SimplifyOptions& simplify_opts
		      , const bool polarity_cnf
		      , const bool notless
		      , const bool input_cuts_only
		      , const bool permute_cnf
		      , const unsigned int permute_cnf_seed
		      )
{
  internal_error("no MiniSAT included");
  exit(1);
}


#else
/*
 *
 * MiniSAT is included
 *
 */
/* Some Minisat include files */
#if defined(MINISAT220CORE)
#include "Solver.h"
#elif defined(MINISAT220SIMP)
#include "SimpSolver.h"
#else
#error "Unknown MiniSAT version defined"
#endif


int BC::minisat_solve(const bool perform_simplifications
		      , const SimplifyOptions& simplify_opts
		      , const bool polarity_cnf
		      , const bool notless
		      , const bool input_cuts_only
		      , const bool permute_cnf
		      , const unsigned int permute_cnf_seed
		      )
{
  bool result;
  int max_var_num;
  unsigned int max_clause_length;
  unsigned int nof_clauses = 0;
#if defined(MINISAT220CORE)
  Minisat::Solver *solver = 0;
#elif defined(MINISAT220SIMP)
  Minisat::SimpSolver *solver = 0;
#else
#error "Unknown MiniSAT version defined"
#endif

  Minisat::Var *map_gatenum_to_minisat_var = 0;

  Timer timer;

  if(perform_simplifications)
    {
      if(!simplify(simplify_opts))
	return 0;
    }
  else
    {
      if(!share())
	return 0;
    }
  

  if(!cnf_normalize())
    return 0;
  
  if(perform_simplifications)
    {
      SimplifyOptions _opts = simplify_opts;
      _opts.preserve_cnf_normalized_form = true;
      if(!simplify(_opts))
	return 0;
    } 
  else
    {
      if(!share())
	return 0;
    }
    


  
  /*
   * Print some statistics
   */
  if(verbose) {
    unsigned int max_min_height, max_max_height;
    compute_stats(max_min_height, max_max_height);
    verbose_print("Preprocessing time: %.2lf\n", timer.get_duration());
    verbose_print("The max-min height of the circuit is %u\n", max_min_height);
    verbose_print("The max-max height of the circuit is %u\n", max_max_height);
  }

  /* Next compute CNF translation time */
  timer.reset();

  /*
   * Find the relevant gates and number them in temp field
   */
  reset_temp_fields(-1);
  int nof_relevant_gates = 0;
  for(Gate *gate = first_gate; gate; gate = gate->next)
    {
      assert(gate->type != Gate::tTRUE ||
	     (gate->determined && gate->value == true));
      assert(gate->type != Gate::tFALSE ||
	     (gate->determined && gate->value == false));
      if(simplify_opts.use_coi == false or
	 (gate->determined and !gate->is_justified()))
	gate->mark_coi(nof_relevant_gates);
    }
  verbose_print("The circuit has %d relevant gates\n", nof_relevant_gates);
  if(nof_relevant_gates == 0)
    {
      goto sat_exit;
    }
  
  
  /*
   * Renumber the gates in temp-fields and compute
   * (i) maximum clause length and
   * (ii) the number of relevant input gates
   */
  {
    max_clause_length = 2;
    int gate_num = 1;
    int nof_relevant_input_gates = 0;
    for(Gate *gate = first_gate; gate; gate = gate->next)
      {
        if(gate->temp == -1) {
          /* Not relevant */
          continue;
        }
        if(notless and gate->type == Gate::tNOT) {
          /* NOT-less translation */
          assert(!gate->determined);
          assert(gate->children->child->type != Gate::tNOT);
          gate->temp = -1;
        }
        else {
          gate->temp = gate_num++;
        }

        /* Compute max clause length */
        if(gate->count_children() + 1 > max_clause_length)
          max_clause_length = gate->count_children() + 1;

        if(gate->type == Gate::tVAR)
          nof_relevant_input_gates++;
      }
    max_var_num = gate_num;

    verbose_print("The circuit has %d relevant input gates\n", nof_relevant_input_gates);
  }

  /*
   * Permute the CNF variable indices if required
   */
  if(permute_cnf)
    {
      verbose_print("Permuting the CNF variables...");

      unsigned int * const perm = my_perm(max_var_num-1, permute_cnf_seed);
      for(Gate *gate = first_gate; gate; gate = gate->next)
	{
	  assert(gate->temp <= max_var_num-1);
	  if(gate->temp <= 0)
	    continue;
	  gate->temp = perm[gate->temp];
	  assert(gate->temp > 0 && gate->temp <= max_var_num-1);
	}
      free(perm);

      verbose_print("done\n");
    }


  /*
   * Init Minisat
   */
#if defined(MINISAT220CORE)
  solver = new Minisat::Solver();
#elif defined(MINISAT220SIMP)
  solver = new Minisat::SimpSolver();
#else
#error "Unknown MiniSAT version defined"
#endif

  /*
   * Get a minisat variable for each relevant gate
   */
  map_gatenum_to_minisat_var = (Minisat::Var*)calloc(max_var_num, sizeof(Minisat::Var));
  for(int i = 1; i < max_var_num; i++)
    {
      map_gatenum_to_minisat_var[i] = solver->newVar();
    }


  /*
   * Compute polarity info if needed
   */
  if(polarity_cnf)
    mir_compute_polarity_information();



  /*
   * Build and feed the CNF to MiniSat
   */
  {
    Minisat::vec<Minisat::Lit> clause;
    std::list<std::vector<int> *> clauses;
    for(Gate *gate = first_gate; gate; gate = gate->next)
      {
	assert(gate->temp == -1 || (gate->temp>0 && gate->temp<max_var_num));
        if(gate->temp == -1) {
          /* Not relevant */
          continue;
        }
        /*
         * Get clauses
         */
	if(polarity_cnf)
	  gate->cnf_get_clauses_polarity(clauses, notless);
	else
	  gate->cnf_get_clauses(clauses, notless);

        /*
         * Add clauses to Minisat
         */
        while(!clauses.empty())
	  {
	    std::vector<int> *cl = clauses.back();
	    clauses.pop_back();
	    assert(cl->size() <= max_clause_length);
	    clause.clear();
	    /* Transform clause into Minisat form */
	    for(std::vector<int>::iterator li = cl->begin();
		li != cl->end();
		li++)
	      {
		int lit = *li;
		assert(lit != 0 && abs(lit) < max_var_num);
		Minisat::Lit minisat_lit = Minisat::mkLit(map_gatenum_to_minisat_var[abs(lit)]);
		if(lit < 0)
		  minisat_lit = ~minisat_lit;
		clause.push(minisat_lit);
	      }
	    /* Add clause to Minisat */
	    solver->addClause(clause);
	    nof_clauses++;
      delete cl;/*AM : ADDED TO AVOID MEMORY LEAKS*/
	  }
	/*
         * Add unit clauses for constrained gates
         */
        if(gate->determined)
	  {
	    clause.clear();
	    Minisat::Lit minisat_lit = Minisat::mkLit(map_gatenum_to_minisat_var[gate->temp]);
	    bool negated = false;
	    if((gate->value == false) ^ negated)
	      minisat_lit = ~minisat_lit;
	    clause.push(minisat_lit);
	    solver->addClause(clause);
	    nof_clauses++;
	  }
	else
	  {
	    /* The following cases should really not happen... */
	    if(gate->type == Gate::tTRUE)
	      {
		clause.clear();
		clause.push(Minisat::mkLit(map_gatenum_to_minisat_var[gate->temp]));
		solver->addClause(clause);
		nof_clauses++;
	      }
	    else if(gate->type == Gate::tFALSE)
	      {
		clause.clear();
		clause.push(~Minisat::mkLit(map_gatenum_to_minisat_var[gate->temp]));
		solver->addClause(clause);
		nof_clauses++;
	      }
	  }
      }
  }

  /*
   * Mark branchable variables
   */
  if(input_cuts_only)
    {
      for(Gate *gate = first_gate; gate; gate = gate->next)
        {
	  assert(gate->temp == -1 or
		 (gate->temp > 0 and gate->temp < max_var_num));
	  /* Not relevant? */
	  if(gate->temp == -1)
	    continue;
	  /* An input or constant gate? */
          if(gate->type == Gate::tVAR or
	     gate->type == Gate::tFALSE or
	     gate->type == Gate::tTRUE)
	    {
	      ;
	    }
	  else
	    {
	      /* Disable branching on this gate */
	      solver->setDecisionVar(map_gatenum_to_minisat_var[gate->temp],
				     false);
	    }
	}
    }

  verbose_print("CNF translation time: %.2lf\n", timer.get_duration());
  verbose_print("The cnf has %d variables and %d clauses\n",
		max_var_num-1, nof_clauses);


  /*
   * Execute minisat
   */
  verbose_print("Executing minisat...\n");

  /* Next measure time spent in Minisat */
  timer.reset();
  solver->verbosity = 2;
  result = solver->solve();
  
  if(verbose) {
    verbose_print("Minisat time: %.2lf\n", timer.get_duration());
    verbose_print("Minisat statistics:\n");
    verbose_print("restarts              : %lu\n",
		  solver->starts);
    verbose_print("conflicts             : %-12lu\n",
		  solver->conflicts);
    verbose_print("decisions             : %-12lu\n",
		  solver->decisions);
    verbose_print("propagations          : %-12lu\n",
		  solver->propagations);
    verbose_print("conflict literals     : %-12lu   (%4.2f %% deleted)\n", solver->tot_literals, (solver->max_literals - solver->tot_literals)*100 / (double)solver->max_literals);
  }

  
  if(result == false)
    {
      free(map_gatenum_to_minisat_var); map_gatenum_to_minisat_var = 0;
      delete solver; solver = 0;
      return 0;
    }

 sat_exit:
  if(solver)
    {
      /*
       * Move Minisat truth assignment on input variables to the circuit
       */
      for(Gate* gate = first_gate; gate; gate = gate->next)
	{
	  assert(gate->temp == -1 or
		 (gate->temp > 0 and gate->temp < max_var_num));
	  if(gate->temp <= 0)
	    continue;
	  if(gate->type != Gate::tVAR)
	    continue;
	  Minisat::lbool val = solver->model[map_gatenum_to_minisat_var[gate->temp]];
	  assert(val == Minisat::lbool(false) or val == Minisat::lbool(true));
	  bool negated = false;
	  const bool minisat_value = (val == Minisat::lbool(true)) ^ negated;
	  if(gate->determined)
	    {
	      if(gate->value != minisat_value)
		{
		  fprintf(stderr, "Inernal error: solution is inconsistent\n");
		  exit(-1);
		}
	    }
	  else
	    {
	      gate->determined = true;
	      gate->value = minisat_value;
	    }
	}
      free(map_gatenum_to_minisat_var); map_gatenum_to_minisat_var = 0;
      delete solver; solver = 0;
    }
  
  /*
   * Assign irrelevant input gates to arbitrary values
   */
  for(Gate* gate = first_gate; gate; gate = gate->next)
    {
      if(gate->type == Gate::tVAR and !gate->determined)
	{
	  assert(gate->temp == -1);
	  gate->determined = true;
	  gate->value = false;
	}
    }
  
  /*
   * Evaluate rest of the irrelevant gates
   */
  for(Gate* gate = first_gate; gate; gate = gate->next)
    {
      if(!gate->determined)
	{
	  const bool evaluation_ok = gate->evaluate();
	  if(!evaluation_ok)
	    internal_error("%s:%u: Evaluation error",__FILE__,__LINE__);
	  DEBUG_ASSERT(gate->determined);
	}
    }

  /*
   * Check consistency
   */
  if(!check_consistency())
    {
      internal_error("%s:%u: Consistency check failed",__FILE__,__LINE__);
      exit(-1);
    }

  return 1;
}


#endif //BC_HAS_MINISAT
