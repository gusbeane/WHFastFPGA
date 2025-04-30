/**
 * @file 	communication_mpi.c
 * @brief	Handles communication between nodes using MPI.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details These routines handle the communication between
 * different nodes via the Message Passing Interface (MPI). 
 * There are two different types of communications implemented 
 * at the moment:
 * - Distributing particles to the correct node.
 * - Creating, and distributing the essential tree to allow
 *   other nodes walk remote trees. Note that the opening 
 *   criteria is different for gravity and collision 
 *   tree walks.
 * 
 * 
 * @section LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
