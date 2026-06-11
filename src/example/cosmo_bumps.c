/* Two-channel H2 formation through cosmic expansion (z=800 -> 10), showing the
   TWO bumps as the CMB stops photo-dissociating H2+ (first, z~400-600) and then
   stops photo-detaching H- (second, z~100-150).  Galli & Palla (1998).

   The recombination + thermal history (x_e(z), T(z)) is IMPOSED from CICASS's
   RECFAST table each step (Grackle is not a recombination code); Grackle's network
   evolves only the molecular intermediaries (H-, H2+, H2, HD) with the CMB
   photo-rates (our fork). This is the correct division of labour and the only way
   the CMB suppression -> release sequence (the two bumps) is resolved.
   Usage: ./cosmo_bumps > cosmo_bumps.dat */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grackle.h>

#define mh 1.67262171e-24
static double H0_s, Om;
static double cosmic_time(double z){ return (2.0/3.0)/H0_s/sqrt(Om)*pow(1.0+z,-1.5); }

/* RECFAST table (z, x_e, T_gas), z descending; interpolate. */
static int NR; static double RZ[300], RX[300], RT[300];
static void load_recfast(const char*p){
  FILE*f=fopen(p,"r"); if(!f){fprintf(stderr,"recfast not found: %s\n",p);exit(1);}
  char line[256]; int i=0; if(!fgets(line,256,f)){exit(1);}  /* count line */
  while(fgets(line,256,f)){ double z,x,t; if(sscanf(line,"%lf %lf %lf",&z,&x,&t)==3){RZ[i]=z;RX[i]=x;RT[i]=t;i++;} }
  fclose(f); NR=i;   /* table is z descending: RZ[0]=1630 ... RZ[NR-1]=0 */
}
static double rf_interp(double z,double*Y){
  if(z>=RZ[0]) return Y[0]; if(z<=RZ[NR-1]) return Y[NR-1];
  for(int i=0;i<NR-1;i++){ if(z<=RZ[i]&&z>=RZ[i+1]){ double w=(z-RZ[i])/(RZ[i+1]-RZ[i]); return Y[i]*(1-w)+Y[i+1]*w; } }
  return Y[NR-1];
}

int main(void){
  double h=0.71; Om=0.27; double Ob=0.046, XH=0.76, DtoH=2.6e-5;
  H0_s=h*100.0*1e5/3.0857e24; double nH0=Ob*1.8788e-29*h*h*XH/mh;
  load_recfast("/Users/tabel/Projects/cicass/vbc_transfer/recfast/xeTrecfast.out");

  code_units u; u.comoving_coordinates=0; u.density_units=mh;
  u.length_units=3.0857e21; u.time_units=3.1557e13; u.a_units=1.0;
  double z=800.0; u.a_value=1.0/(1.0+z); set_velocity_units(&u);

  chemistry_data*cd=malloc(sizeof(chemistry_data)); set_default_chemistry_parameters(cd);
  grackle_data->use_grackle=1; grackle_data->with_radiative_cooling=1;
  grackle_data->primordial_chemistry=3; grackle_data->metal_cooling=0;
  grackle_data->dust_chemistry=0; grackle_data->UVbackground=0;
  grackle_data->cmb_dissociation=1;
  grackle_data->grackle_data_file="../../input/CloudyData_noUVB.h5";
  if(initialize_chemistry_data(&u)==0){fprintf(stderr,"init failed\n");return 1;}

  grackle_field_data f; int dim[3]={1,1,1},st[3]={0,0,0},en[3]={0,0,0};
  f.grid_rank=3; f.grid_dimension=dim; f.grid_start=st; f.grid_end=en; f.grid_dx=0;
  gr_float dens[1],ie[1],vx[1]={0},vy[1]={0},vz[1]={0};
  gr_float HI[1],HII[1],HeI[1],HeII[1],HeIII[1],e[1],HM[1],H2I[1],H2II[1],DI[1],DII[1],HDI[1],Z[1];
  f.density=dens; f.internal_energy=ie; f.x_velocity=vx; f.y_velocity=vy; f.z_velocity=vz;
  f.HI_density=HI; f.HII_density=HII; f.HeI_density=HeI; f.HeII_density=HeII; f.HeIII_density=HeIII;
  f.e_density=e; f.HM_density=HM; f.H2I_density=H2I; f.H2II_density=H2II;
  f.DI_density=DI; f.DII_density=DII; f.HDI_density=HDI; f.metal_density=Z;
  f.volumetric_heating_rate=NULL; f.specific_heating_rate=NULL;
  f.RT_HI_ionization_rate=NULL; f.RT_HeI_ionization_rate=NULL; f.RT_HeII_ionization_rate=NULL;
  f.RT_H2_dissociation_rate=NULL; f.RT_heating_rate=NULL;

  double tiny=1e-20, Tunits=get_temperature_units(&u), mu=1.22;
  double rho=nH0*pow(1+z,3.0)*mh/XH; dens[0]=rho/u.density_units;
  HM[0]=tiny*dens[0]; H2I[0]=tiny*dens[0]; H2II[0]=tiny*dens[0];
  DII[0]=tiny*dens[0]; HDI[0]=tiny*dens[0]; Z[0]=tiny*dens[0];

  printf("# z T_gas T_cmb xe xH2 xHM xH2p xHD\n");
  double dlna=0.002;
  while(z>10.0){
    double zn=(1.0+z)*exp(-dlna)-1.0; if(zn<10.0) zn=10.0;
    double f3=pow((1.0+zn)/(1.0+z),3.0);
    dens[0]*=f3; HM[0]*=f3; H2I[0]*=f3; H2II[0]*=f3; DII[0]*=f3; HDI[0]*=f3; Z[0]*=f3;
    /* IMPOSE recfast ionization + thermal history; keep the molecules and conserve
       H nuclei (HI = total H - H+ - the H locked in H-, H2, H2+). */
    double xe=rf_interp(zn,RX), Tg=rf_interp(zn,RT);
    HII[0]=xe*XH*dens[0]; e[0]=xe*XH*dens[0];
    HI[0]=XH*dens[0] - HII[0] - HM[0] - H2I[0] - H2II[0];
    if(HI[0]<tiny*dens[0]) HI[0]=tiny*dens[0];
    HeI[0]=(1-XH)*dens[0]; HeII[0]=tiny*dens[0]; HeIII[0]=tiny*dens[0];
    DI[0]=(1-xe)*DtoH*XH*dens[0];
    ie[0]=Tg/Tunits/(5.0/3.0-1.0)/mu;
    u.a_value=1.0/(1.0+zn);
    double dt=(cosmic_time(zn)-cosmic_time(z))/u.time_units;
    if(solve_chemistry(&u,&f,dt)==0){fprintf(stderr,"solve failed z=%g\n",z);return 1;}
    double nHn=HI[0]+HII[0]+HM[0]+H2I[0]+H2II[0];
    printf("%.4f %.3e %.3e %.4e %.4e %.4e %.4e %.4e\n", zn, Tg, 2.73*(1+zn),
           xe, (H2I[0]/2.0)/nHn, HM[0]/nHn, (H2II[0]/2.0)/nHn, (HDI[0]/3.0)/nHn);
    z=zn;
  }
  return 0;
}
