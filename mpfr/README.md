Vendored dependency: MPFR C++ (mpreal)
======================================

    File:              mpreal.h
    Version:           3.7.2
    Upstream:          https://github.com/advanpix/mpreal
    Upstream commit:   b864ace57f53  (2025-11-22)
    Imported:          2026-08-21
    Local changes:     none — byte-for-byte upstream

Update procedure
----------------

    git clone https://github.com/advanpix/mpreal /tmp/mpreal
    cp /tmp/mpreal/mpreal.h mpfr/mpreal.h
    # then update the four fields above

Before importing, confirm the current file is unmodified, so a local fix is not
silently discarded:

    diff <(tr -d '\r' < mpfr/mpreal.h) <(git -C /tmp/mpreal show <recorded-sha>:mpreal.h)

Why this record exists
----------------------

This header sat at 3.6.2 (upstream commit 0370b406, 24 July 2015) for a decade
with nothing recording where it came from or when. It stopped compiling against
MPFR 4.x, which is what current distributions ship, and the first CI run in ten
years died on it:

    mpfr/mpreal.h:582:27: error: no matching function for call to
        'mpfr::mpreal::mpfr_srcptr(mpfr_srcptr)'

Version 3.7.2 carries explicit MPFR_VERSION_NUM(4,0,0) guards. Since
HJCFIT_USE_MPFR defaults to ON for UNIX, this header is in the default build
path, not an optional extra — MPFR provides the arbitrary-precision fallback for
root bracketing in likelihood/root_finder.cc when double precision fails.

Only mpfr::mpreal and mpreal::set_default_prec are used directly; Eigen
integration comes from Eigen's own unsupported/Eigen/MPRealSupport, included
via HJCFITConfig.h.in.

--------------------------------------------------------------------------

Upstream description and licence
================================

    MPFR C++ (MPREAL):     Multiple precision floating point arithmetic library for C++.
                           Thread-safe, cross-platform (MSVC, GCC, ICC), one-header C++ library.
                           Supports C++11 features if available, C++03 compatible otherwise.
                           
    Thin wrapper for MPFR: http://mpfr.org
	
    Project homepage:    http://www.holoborodko.com/pavel/mpfr
    Contact e-mail:      pavel@holoborodko.com
    
    Copyright (c) 2008-2015 Pavel Holoborodko

    Contributors:
    Dmitriy Gubanov, Konstantin Holoborodko, Brian Gladman, 
    Helmut Jarausch, Fokko Beekhof, Ulrich Mutze, Heinz van Saanen, 
    Pere Constans, Peter van Hoof, Gael Guennebaud, Tsai Chia Cheng, 
    Alexei Zubanov, Jauhien Piatlicki, Victor Berger, John Westwood,
    Petr Aleksandrov, Orion Poplawski, Charles Karney, Arash Partow,
    Rodney James, Jorge Leitao.

    Licensing:
    (A) MPFR C++ is under GNU General Public License ("GPL").
    
    (B) Non-free licenses may also be purchased from the author, for users who 
        do not want their programs protected by the GPL.

        The non-free licenses are for users that wish to use MPFR C++ in 
        their products but are unwilling to release their software 
        under the GPL (which would require them to release source code 
        and allow free redistribution).

        Such users can purchase an unlimited-use license from the author.
        Contact us for more details.
        
    GNU General Public License ("GPL") copyright permissions statement:
    **************************************************************************
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
