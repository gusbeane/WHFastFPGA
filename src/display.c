/**
 * @file    display.c
 * @brief   Realtime OpenGL visualization.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details These functions provide real time visualizations
 * using OpenGL. 
 * 
 * @section LICENSE
 * Copyright (c) 2015 Hanno Rein, Shangfei Liu
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
#define DEG2RAD (M_PI/180.)
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include "rebound.h"
#include "display.h"
#include "tools.h"
#include "particle.h"
#include "boundary.h"
#include "display.h"
#include "output.h"
#include "integrator.h"
#define MAX(a, b) ((a) < (b) ? (b) : (a))       ///< Returns the maximum of a and b

static void reb_display_set_default_view(struct reb_simulation* const r, struct reb_display_settings* s){
    float  scale = 0.;
    // Need a scale for visualization
    if (r->root_size==-1){  
        scale = 0.;
        const struct reb_particle* p = r->particles;
        for (unsigned int i=0;i<r->N-r->N_var;i++){
            const double _r = sqrt(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z);
            scale = MAX(scale, _r);
        }
        if(scale==0.){
            scale = 1.;
        }
        scale *= 1.1;
    }else{
        scale = r->boxsize_max/2.;
    }
    
    struct reb_mat4df oldview = s->view;
    s->view = reb_mat4df_scale(reb_mat4df_identity(), 1./scale, 1./scale, 1./scale);
    if (oldview.m[1]==0. && oldview.m[2]==0. && oldview.m[4]==0. && oldview.m[6]==0.){
        struct reb_rotation rotation = {
            .ix = 1./sqrt(2.),
            .iy = 0.,
            .iz = 0.,
            .r = 1./sqrt(2.),
        };
        s->view = reb_mat4df_multiply(reb_rotation_to_mat4df(rotation), s->view);
    }else if (oldview.m[1]==0. && oldview.m[2]==0. && oldview.m[4]==0. && oldview.m[5]==0.){
        struct reb_rotation rotation = {
            .ix = 0.,
            .iy = -1./sqrt(2.),
            .iz = 0.,
            .r = 1./sqrt(2.),
        };
        s->view = reb_mat4df_multiply(reb_rotation_to_mat4df(rotation), s->view);
    }
}

void reb_display_settings_init(struct reb_simulation*r, struct reb_display_settings* s){
    if (r->max_radius0 > 0.0){
        s->spheres       = 1; 
    }else{
        s->spheres       = 0; 
    }
    s->pause             = 0; 
    s->multisample       = 1; 
    if (r->integrator==REB_INTEGRATOR_WHFAST){
        s->wire          = 1; 
    }else{
        s->wire          = 0; 
    }
    s->breadcrumbs       = 0;
    s->onscreentext      = 1; 
    s->ghostboxes        = 0; 
    s->reference         = -1;
    s->view.m[1]=1; // this will make set_default_view show the xy plane
    reb_display_set_default_view(r, s);
}

void reb_simulation_add_display_settings(struct reb_simulation*r){
    if (r->display_settings){
        reb_simulation_error(r,"Simulation already has display settings.");
        return;
    }
    r->display_settings = calloc(1,sizeof(struct reb_display_settings));
    reb_display_settings_init(r, r->display_settings);
}




