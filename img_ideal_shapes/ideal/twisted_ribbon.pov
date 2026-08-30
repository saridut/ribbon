#version 3.7;
global_settings{assumed_gamma 1.0}

#include "colors.inc"
#include "textures.inc"
#include "metals.inc"
#include "glass.inc"
#include "rad_def.inc"
#include "twisted_ribbon.inc" 

#default{finish{ambient 0.5 diffuse 0.9}}


global_settings{
	radiosity{
		Rad_Settings(Radiosity_IndoorHQ, on, on)
	}	
}

global_settings {
  photons {
    count 20000
    autostop 0
    jitter .4
    }
  }
  


object {
TwistedRibbon    
texture { T_Brass_1C }
photons{
  target
  reflection on
  refraction on
  }
rotate <0, 90, 0>
rotate <-15, 0, 0>
rotate <0, -30, 0>
}


// the floor along the x-z plane (y is the normal vector)
plane { 
  y, -15  
  pigment { color White }   
}


background {color White*0.3}

light_source { <0, 25, -5> color White 
	area_light <5, 0, 0>, <0, 0, 3>, 15, 15
    adaptive 2 jitter}

camera {
  location <0, 18, -13>
  look_at <0, 0, -3>
}

