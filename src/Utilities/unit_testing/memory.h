/*
*************** License ***************
thunder-tester is a simple system for in code unit testsing in C using gcc
Copyright (C) 2019  Tor Djärv email: tordjarv@gmail.com
This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published b
    the Free Software Foundation, either version 3 of the License, o
    (at your option) any later version

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty o
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
GNU General Public License for more details
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef __MEMORY__
#define __MEMORY__
#define MALLOC(size,type) (type*)malloc(size*sizeof(type))
#define CALLOC(size,type) (type*)calloc(size,sizeof(type))
#define REALLOC(array,size) (typeof(array))realloc(array,size*sizeof(typeof(*array)))
#endif
