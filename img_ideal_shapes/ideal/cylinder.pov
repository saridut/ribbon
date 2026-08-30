#version 3.7;
global_settings{assumed_gamma 1.0}

#include "colors.inc"
#include "textures.inc"
#include "metals.inc"

#default{finish{ambient 0.5 diffuse 0.9}}

#include "rad_def.inc"
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
  
  
#declare CsgDelta = 0.01;
#declare Rin = 1.8;
#declare Thickness = 0.1;
#declare Gap = 0.4;
#declare Rout = Rin + Thickness;
#declare Length = 20.0;

difference {
cylinder {
    < 0, 0, -0.5*Length >, < 0, 0, 0.5*Length >, Rout
}
cylinder {
    < 0, 0, -(0.5*Length+CsgDelta) >, < 0, 0, (0.5*Length+CsgDelta) >, (Rout-Thickness)
}
box{
	<-0.5*Gap, Rout-2*Thickness, -(0.5*Length+2*CsgDelta) >, 
	< 0.5*Gap, Rout+2*Thickness,  (0.5*Length+2*CsgDelta) >
}

texture { T_Brass_1C }
photons{
  target
  reflection on
  refraction on
  }

translate 0.5*y
rotate <0, 105, 0>
rotate <-10, 0, 0>
}


// the floor along the x-z plane (y is the normal vector)
plane { 
  y, -15  
  pigment { color White }   
}


background {color White*0.3}

light_source { <3, 25, -4> color White 
	area_light <5, 0, 0>, <0, 0, 3>, 15, 15
    adaptive 2 jitter}

camera {
  location <0, 20, -18>
  look_at <0, 0, -3>
}

