
// Just the code for the color bar.

/* List of choices:
0. Original
1. Standard
2. Copper
3. Roy G. Biv
4. Electric
5. Blank
6. Inv. Blank
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void get_rgb( double val , float * rp , float * gp , float * bp , int COLORBAR ){
   
   double rrr,ggg,bbb;
   double v5, v4, v3, v2;
   v2 = val*val;
   v3 = val*v2;
   v4 = v2*v2;
   v5 = v3*v2;

   if( COLORBAR == 0 ){
      double nexp = 8.0;
      rrr = exp(-nexp*pow(val-5./6.,2.0)) + .25*exp(-nexp*pow(val+1./6.,2.0));
      ggg = exp(-nexp*pow(val-3./6.,2.0));
      bbb = exp(-nexp*pow(val-1./6.,2.0)) + .25*exp(-nexp*pow(val-7./6.,2.0));
   }else if(COLORBAR == 1){
      if( val < .1 ){
         bbb = 4.*(val+.15);
         ggg = 0.0;
         rrr = 0.0;
      }else if( val < .35){
         bbb = 1.0;
         ggg = 4.*(val-.1);
         rrr = 0.0;
      }else if( val < .6 ){
         bbb = 4.*(.6-val);
         ggg = 1.;
         rrr = 4.*(val-.35);
      }else if( val < .85){
         bbb = 0.0;
         ggg = 4.*(.85-val);
         rrr = 1.;
      }else{
         bbb = 0.0;
         ggg = 0.0;
         rrr = 4.*(1.1-val);
      }
   }else if(COLORBAR == 2){
      rrr = 2.*val;
      ggg = 1.2*val;
      bbb = .8*val;
   }else if(COLORBAR == 3){
      double gam = .8;
      double Amp;
      double r0,g0,b0;
      double hi,lo,x1,x2,x3,x4;
      hi = .8;
      lo = .1;
      if( val > hi ) Amp = .3 + .7*(1.-val)/(1.-hi);
      else if( val < lo ) Amp = .3 + .7*(val)/(lo);
      else Amp = 1.0;

      x1 = .5;
      x2 = .325;
      x3 = .15;
      x4 = 0.;

      if( val > x1 )      r0 = 1.;
      else if( val > x2 ) r0 = (val-x2)/(x1-x2);
      else if( val > x3 ) r0 = 0.;
      else if( val > x4 ) r0 = (val-x3)/(x4-x3);
      else                r0 = 1.;

      x1 = .6625;
      x2 = .5;
      x3 = .275;
      x4 = .15;

      if( val > x1 )      g0 = 0.;
      else if( val > x2 ) g0 = (val-x1)/(x2-x1);
      else if( val > x3 ) g0 = 1.;
      else if( val > x4 ) g0 = (val-x4)/(x3-x4);
      else                g0 = 0.;

      x1 = .325;
      x2 = .275;

      if( val > x1 )      b0 = 0.;
      else if( val > x2 ) b0 = (val-x1)/(x2-x1);
      else                b0 = 1.;

      rrr = pow(Amp*r0,gam);
      ggg = pow(Amp*g0,gam);
      bbb = pow(Amp*b0,gam);
      
      //if( val > .99 ){ rrr = 1.0; ggg = 1.0 ; bbb = 1.0; }
      //if( val < .01 ){ rrr = 0.0; ggg = 0.0 ; bbb = 0.0; }

   }else if(COLORBAR == 4){
      if( val < .1 ){
         bbb = 4.*(val+.125);
         ggg = 0.0;
         rrr = 0.0;
      }else if( val < .375){
         bbb = 1.0;
         ggg = 4.*(val-.125);
         rrr = 0.0;
      }else if( val < .625 ){
         bbb = 4.*(.625-val);
         rrr = 4.*(val-.375);
         ggg = bbb;
         if( rrr > bbb ) ggg = rrr;
      }else if( val < .875){
         bbb = 0.0;
         ggg = 4.*(.875-val);
         rrr = 1.;
      }else{
         bbb = 0.0;
         ggg = 0.0;
         rrr = 4.*(1.125-val);
      }
   }else if(COLORBAR == 5){
      rrr = val;
      ggg = val;
      bbb = val;
   }else if(COLORBAR == 6){
      rrr = 1.-val;
      ggg = 1.-val;
      bbb = 1.-val;
   }else if(COLORBAR == 7){
      // magma
      rrr = 5.2687848*v5 - 11.44778744*v4  + 6.13567269*v3 -0.01964937*v2 + 1.0698287*val -0.02223323;
      ggg = -5.54509553*v5 + 11.4659551*v4 - 6.76338190*v3 + 1.63850072*v2 + 0.178140470*val + .0105032349;
      bbb = -12.8819794*v5 + 32.9740674*v4 -24.7317822*v3 + 3.19023254*v2 + 2.21260079*val -0.00661976479;
      if (rrr < 0.0) rrr = 0.0;
      if (rrr > 1.0) rrr = 1.0;
      if (ggg < 0.0) ggg = 0.0;
      if (ggg > 1.0) ggg = 1.0;
      if (bbb < 0.0) bbb = 0.0;
      if (bbb > 1.0) bbb = 1.0;
   }else if(COLORBAR == 8){
      // inferno
      double v6, v7, v8, v9;
      v6 = v3*v3;
      v7 = v3*v4;
      v8 = v4*v4;
      v9 = v4*v5;
      rrr = 353.03327105*v9 + -1470.88275954*v8 + 2518.04284713*v7 + -2253.31553816*v6 + 1093.49189644*v5 + -254.27368291*v4 + 5.51981471*v3 + 9.33740110*v2 + 0.03763654*val + 0.00228385;
      ggg = -558.70563289*v9 + 2429.27311490*v8 + -4344.49539941*v7 + 4094.43416755*v6 + -2155.38216084*v5 + 609.75414705*v4 + -74.45594078*v3 + -0.22560533*v2 + 0.80290083*val + -0.00815033;
      bbb = 434.51062082*v9 + -1934.52177132*v8 + 3379.24756989*v7 + -2881.66653015*v6 + 1156.34385196*v5 + -99.43855605*v4 + -70.96109647*v3 + 15.80548478*v2 + 1.31196143*val + 0.01188670;
      if (rrr < 0.0) rrr = 0.0;
      if (rrr > 1.0) rrr = 1.0;
      if (ggg < 0.0) ggg = 0.0;
      if (ggg > 1.0) ggg = 1.0;
      if (bbb < 0.0) bbb = 0.0;
      if (bbb > 1.0) bbb = 1.0;
   }else if(COLORBAR == 9){
      // viridis
      rrr = -11.80729026*v5 + 25.27758373*v4 -14.85075402*v3  + 2.2290938*v2 -0.14187082*val + 0.2800557;
      ggg = 0.121926730*v5 -1.52661275*v4 + 2.53827574*v3 -1.83476402*v2 + 1.60267540*val -0.00126998147;
      bbb = 13.69761657*v5 -33.28507238*v4 +  28.69652089*v3 -11.89994842*v2 + 2.58966906*val + 0.30238385;
      if (rrr < 0.0) rrr = 0.0;
      if (rrr > 1.0) rrr = 1.0;
      if (ggg < 0.0) ggg = 0.0;
      if (ggg > 1.0) ggg = 1.0;
      if (bbb < 0.0) bbb = 0.0;
      if (bbb > 1.0) bbb = 1.0;

   }else if(COLORBAR == 10){
      // fearless idea
      double v6, v7, v8, v9, v10, v11, v12, v13, v14;
      v6 = v3*v3;
      v7 = v3*v4;
      v8 = v4*v4;
      v9 = v4*v5;
      v10 = v5*v5;
      v11 = v5*v6;
      v12 = v6*v6;
      v13 = v6*v7;
      v14 = v7*v7;

      rrr = -339052.58765971*v14 + 2378322.62792338*v13 + -7431459.27726085*v12 + 13634285.65033955*v11 + -16299573.61354742*v10 + 13323480.04072026*v9 + -7603567.35982416*v8 + 3039969.93747006*v7 + -842381.67114325*v6 + 157760.12200660*v5 + -19136.60729466*v4 + 1408.55471235*v3 + -58.53606105*v2 + 3.68703562*val + 0.02679899;
      ggg = 345834.18942290*v14 + -2407536.68231237*v13 + 7466528.29046747*v12 + -13600398.02064648*v11 + 16151959.61619647*v10 + -13128280.10495173*v9 + 7460077.94108592*v8 + -2975253.51067690*v7 + 824368.35915176*v6 + -154913.79382242*v5 + 18987.96057112*v4 + -1432.63857532*v3 + 59.88173696*v2 + -0.49050263*val + 0.00621574;
      bbb = 977692.24248626*v14 + -6885588.27955981*v13 + 21614867.99947388*v12 + -39865230.47268245*v11 + 47936542.22890460*v10 + -39427057.65051076*v9 + 22638900.00981719*v8 + -9099458.94097637*v7 + 2529745.39770873*v6 + -473416.11286858*v5 + 56927.12796437*v4 + -4066.73266942*v3 + 144.67519126*v2 + -0.52597951*val + 0.02192712;

      if (rrr < 0.0) rrr = 0.0;
      if (rrr > 1.0) rrr = 1.0;
      if (ggg < 0.0) ggg = 0.0;
      if (ggg > 1.0) ggg = 1.0;
      if (bbb < 0.0) bbb = 0.0;
      if (bbb > 1.0) bbb = 1.0;
   }else if(COLORBAR == 11){
      // Red-Blue
      rrr = 10.20752512*v5 -17.16433678*v4 + 8.88247458*v3 -5.20403454*v2 + 2.94959219*val + 0.41422865;
      ggg = -11.5674236*v5 + 40.4428654*v4 -49.2169767*v3 +21.4518256*v2 -.941994363*val + 0.0307091725;
      bbb = -8.89193203*v5 + 32.18016432*v4 -42.5403857*v3 + 21.78773187*v2 - 2.3075258*val + 0.19384687;
      if (rrr < 0.0) rrr = 0.0;
      if (rrr > 1.0) rrr = 1.0;
      if (ggg < 0.0) ggg = 0.0;
      if (ggg > 1.0) ggg = 1.0;
      if (bbb < 0.0) bbb = 0.0;
      if (bbb > 1.0) bbb = 1.0;
    }else if(COLORBAR == 12){
      // Purple-Green
      double v6, v7;
      v6 = v3*v3;
      v7 = v3*v4;
      rrr = -3.15151489*v7 + -26.44357655*v6 + 111.55029200*v5 + -141.32019992*v4 + 72.90100974*v3 + -17.05941751*v2 + 3.29985392*val + 0.23189982;
      ggg = -37.98507357*v7 + 129.97007419*v6 + -174.03829563*v5 + 121.99570095*v4 + -53.84964418*v3 + 13.33984287*v2 + 0.82045038*val + 0.00567647;
      bbb = -51.85314990*v7 + 134.68646793*v6 + -103.62306658*v5 + 6.59218715*v4 + 19.59112571*v3 + -8.49746250*v2 + 2.91776116*val + 0.27960782;
      if (rrr < 0.0) rrr = 0.0;
      if (rrr > 1.0) rrr = 1.0;
      if (ggg < 0.0) ggg = 0.0;
      if (ggg > 1.0) ggg = 1.0;
      if (bbb < 0.0) bbb = 0.0;
      if (bbb > 1.0) bbb = 1.0;
  }else{
      rrr = 1.0;
      ggg = 1.0;
      bbb = 1.0;
   }

   *rp = rrr;
   *gp = ggg;
   *bp = bbb;

}
