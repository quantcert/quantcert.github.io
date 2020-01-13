#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module is an simple optimization module that helps find the point of 
maximum value for a given function.
"""
from sage.all import *

def optimize(func, args_init, step_init, step_min, iter_max, verbose=False):
  r"""
  Optimization function finding the maximum reaching coordinates for ``func``
  with a random walk. 
  
  :param function func: The function to be optimized, each of its arguments must 
    be numerical and will be tweaked to find ``func``'s maximum. 
  :param tuple args_init: Initial coordinates for the random walk. 
  :param real step_init: The step size will vary with time in this function, so
    this is the initial value for the step size. 
  :param real step_min: The limit size for the step.
  :param int iter_max: Upper iterations bound for each loop to avoid infinite 
    loops.
  :param bool verbose: If ``verbose`` then extra run information will be 
    displayed in terminal.
  :returns: vector, complex -- The coordinates of the optimum found for 
    ``func`` and the value of ``func`` at this point.
  """
  solution = vector(CC, args_init).normalized()
  solution_temp = solution
  value = func(solution)
  value_temp = value

  step_current = step_init

  if verbose:
    print("Initial solution : " + func.__name__ + 
      str(solution) + " = " + str(value))

  iter_nb = 0
  while step_current > step_min:
    count = 1
    while value_temp <= value and count < iter_max:
      direction = vector(
        [uniform(-1,1) for _ in range(len(args_init))]
        ).normalized()
      solution_temp = solution + direction*step_current
      value_temp = func(solution_temp)
      count += 1
    if value_temp > value:
      solution = solution_temp
      value = value_temp
    else:
      step_current = step_current/2
  
  if verbose:
    print("Final solution : " + func.__name__ + 
      str(solution) + " = " + str(value))

  return (solution, value)


def optimize_normalized(func, normalizer_func, args_init, step_init, step_min, iter_max, verbose=False):
  r"""
  Optimization function finding the maximum reaching coordinates for ``func``
  with a random walk. 
  
  :param function func: The function to be optimized, each of its arguments must 
    be numerical and will be tweaked to find ``func``'s maximum. 
  :param tuple args_init: Initial coordinates for the random walk. 
  :param real step_init: The step size will vary with time in this function, so
    this is the initial value for the step size. 
  :param real step_min: The limit size for the step.
  :param int iter_max: Upper iterations bound for each loop to avoid infinite 
    loops.
  :param bool verbose: If ``verbose`` then extra run information will be 
    displayed in terminal. 
  :returns: vector, complex -- The coordinates of the optimum found for 
    ``func`` and the value of ``func`` at this point.
  """
  solution = vector(CC, normalizer_func(args_init))
  solution_temp = solution
  value = func(solution)
  value_temp = value

  step_current = step_init

  if verbose:
    print("Initial solution : " + func.__name__ + 
      str(solution) + " = " + str(value))

  iter_nb = 0
  while step_current > step_min:
    count = 1
    while value_temp <= value and count < iter_max:
      direction = vector(
        [uniform(-1,1) for _ in range(len(args_init))]
        ).normalized()
      solution_temp = vector(normalizer_func(solution + direction*step_current))
      value_temp = func(solution_temp)
      count += 1
    if value_temp > value:
      solution = solution_temp
      value = value_temp
    else:
      step_current = step_current/2
      if verbose:
        print("optimization, current step: " + str(step_current))
  
  if verbose:
    print("Final solution : " + func.__name__ + 
      str(solution) + " = " + str(value))

  return (solution, value)


def optimize_2spheres(func, args_init, step_init, step_min, iter_max, radius=1, verbose = False):
  r"""
  Optimization function finding the maximum reaching coordinates for ``func``
  with a random walk on a two sphere of dimension half the size of 
  ``args_init``. (Work in progress !)

    For now, this function is in project and is not used, it can be ignored.
  
  :param function func: The function to be optimized, each of its arguments must be 
    numerical and will be tweaked to find ``func``'s maximum. 
  :param tuple args_init: Initial coordinates for the random walk.
  :param real step_init: The step size will vary with time in this function, so
    this is the initial value for the step size. 
  :param real step_min: The limit size for the step.
  :param int iter_max: Upper iterations bound for each loop to avoid infinite 
    loops.
  :param real radius: Sphere radius.
  :param bool verbose: If ``verbose`` then extra run information will be displayed in 
    terminal.
  :returns: vector, complex -- The coordianates of the optimum found for 
    ``func`` and the value of ``func`` at this point.
  """
  def unwrap(vector_tuple):
    """Takes in vectors and returns a list of their coefficients
    """
    arguments_unwraped = []
    for vector_instance in vector_tuple:
      for coefficient in vector_instance:
        arguments_unwraped.append(coefficient)
    return arguments_unwraped

  def point_on_cone(cone_center, cone_spherical_radius, sphere_radius): # TODO
    dimension = len(cone_center)
    rotation_cone_center_to_Z = rotation_to_Z(cone_center) # matrix
    point_random_angles = [uniform(-pi,pi) for _ in dimension - 2]
    cone_angle = cone_spherical_radius/sphere_radius
    point = vector_from_pherical(sphere_radius, cone_angle, *point_random_angles)
    point = rotation_cone_center_to_Z.inverse() * point
    return point

  dimension = len(args_init)/2
  solution = vector(CC, args_init[:dimension]).normalized(), \
    vector(CC, args_init[dimension:]).normalized()
  value = func(unwrap(solution))
  value_temp = value

  step_current = step_init

  if verbose:
    print("Initial solution : " + func.__name__ + 
      str(solution) + " = " + str(value))

  iter_nb = 0
  while step_current > step_min:
    solution_temp = solution # initialization
    count = 1
    while CC(value_temp) <= CC(value) and count < iter_max:
      solution_temp = point_on_cone(solution, step_current, 1)
      value_temp = func(unwrap(solution_temp))
      count += 1
    if CC(value_temp) > CC(value):
      solution = solution_temp
      value = value_temp
    else:
      step_current = step_current/2
  
  if verbose:
    print("Final solution : " + func.__name__ + 
      str(solution) + " = " + str(value))

  return (solution, value)
