//Persistence of Vision Ray Tracer Scene Description File
#version 3.7;
global_settings{assumed_gamma 1.0}

#include "colors.inc"
#include "textures.inc"
#include "rad_def.inc"

#default{finish{ambient 0.5 diffuse 0.9}}

//global_settings{
//	radiosity{
//		Rad_Settings(Radiosity_IndoorLQ, off, off)
//	}	
//}

//Spherical atoms
#declare AtmCd = sphere{<0,0,0>, 1.51 pigment{rgb <204,204,204>/255}}
#declare AtmSe = sphere{<0,0,0>, 1.2  pigment{rgb <255,153,0>  /255}  }
//
#declare FnConfig = "NPL_ML3_105x10_45_caC8.povdat";
#declare NumAtoms = 0;

#fopen FhConfig FnConfig read

#read (FhConfig, NumAtoms)
#debug concat(str(NumAtoms, 10, 0), "\n") 


#declare AtmPos = array[NumAtoms][3]
#declare AtmName = array[NumAtoms]
#declare i = 0;
#for (i, 0, NumAtoms-1)
	#declare AtmPos[i][0] = 0;
	#declare AtmPos[i][1] = 0;
	#declare AtmPos[i][2] = 0;
	#declare AtmName[i] = "00";
#end


//Read atom names and positions
#declare i = 0;
#while (i < NumAtoms)
    #read (FhConfig, AtmName[i], AtmPos[i][0], AtmPos[i][1], AtmPos[i][2])
    #declare i = i + 1;
#end

//#for (i, 0, NumAtoms-1)
//    #read (FhConfig, AtmName[i], AtmPos[i][0], AtmPos[i][1], AtmPos[i][2])
//#end

#fclose FhConfig

//Scene objects

union{
    #for (i, 0, NumAtoms-1)
    	#declare pos = <AtmPos[i][0], AtmPos[i][1], AtmPos[i][2]>;
        #if (AtmName[i] = "Cd")
        	object {AtmCd translate pos}
        #elseif (AtmName[i] = "Se")
        	object {AtmSe translate pos}
        #end
        //finish {phong 1 brilliance 3}
    #end
    rotate <0.0, 0, 0.0>
    scale 1
}

//Background
//background {color White*0.5}

//Lights
light_source { 
    <0, 250, -300> color White 
    area_light <5, 0, 0>, <0, 0, 3>, 15, 15
    adaptive 2 jitter}

//Camera
camera {
  location +10*y-500*z
  look_at  <0, 0, 0> //Put at the end of all camera options
}
