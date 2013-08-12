.. _cpp_occupancies_api:

Occupancies
===========

.. cpp:function:: t_initvec occupancies(MissedEventsG const& _this, bool _initial=true) 

   Solves the equilibrium occupancy equation for initial and final states

   The equilibrium equation is :math:`\phi = \phi M`, :math:`\sum_i [\phi]_i = 1`, where :math:`M` is
   for initial states :math:`\mathcal{G}_{AF}(s=0) \mathcal{G}_{FA}(s=0)`. The problem is solved
   using Eigen's linear least-square utility, adding an extra row to the matrix to impose the second
   condition.

   :param _this: 
      An :cpp:class:`IdealG` or :cpp:class:`MissedEventsG` instance. 
   :param bool _initial:
      Indicates whether to return initial (default, ``true``) or final occupancies.

.. cpp:function:: t_initvec CHS_occupancies(MissedEventsG const& _this, bool _initial=true) 

    Solves the CHS occupancy equation for initial and final states

    For ``_initial=true`` this is Eq 5.11 of :cite:`colquhoun:1996`, whereas for
    ``_initial=false``, it is Eq 5.8.  
    
    .. note::
    
       The initial and final vectors are not quite as satisfyingly symmetric as for other
       occupancies. However, to simplify the API, this function still deals with both cases.

   :param _this: 
      An :cpp:class:`IdealG` or :cpp:class:`MissedEventsG` instance. 
   :param bool _initial:
      Indicates whether to return initial (default, ``true``) or final occupancies.
