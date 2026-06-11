/* Hydrogen recombination history: Grackle (with the v2026 recfast-matched
   recombination: PPB91 alpha_B x Peebles C-factor) vs the RECFAST table.
   One zone, homogeneous, evolved through expansion from z=1000 down.
   CMBREC=1 (default) turns the C-factor on; CMBREC=0 = stock Grackle (which
   over-recombines).  Prints z, x_e(grackle), x_e(recfast), ratio.
   Build: see Makefile target; run from src/example. */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grackle.h>

#define mh 1.67262171e-24
#define kb 1.3806504e-16

/* EdS-ish cosmic time [s] (matter+radiation negligible diff at z<1000 for dt). */
static double cosmic_time(double z, double H0_s, double Om) {
  return (2.0/3.0)/H0_s/sqrt(Om) * pow(1.0+z, -1.5);
}

/* RECFAST table (z descending). */
static int    NT=0;
static double zt[400], xt[400], Tt[400];
static void load_recfast(const char *path) {
  FILE *fp = fopen(path,"r"); if(!fp){fprintf(stderr,"no recfast table %s\n",path); exit(1);}
  char line[256]; fgets(line,sizeof line,fp);           /* count header */
  while (fgets(line,sizeof line,fp)) {
    if (sscanf(line,"%lf %lf %lf",&zt[NT],&xt[NT],&Tt[NT])==3) NT++;
  }
  fclose(fp);
}
static double interp_rf(double z, double *col) {        /* table is z-descending */
  if (z>=zt[0]) return col[0];
  if (z<=zt[NT-1]) return col[NT-1];
  for (int i=0;i<NT-1;i++) if (z<=zt[i] && z>zt[i+1]) {
    double f=(z-zt[i])/(zt[i+1]-zt[i]); return col[i]+f*(col[i+1]-col[i]);
  }
  return col[NT-1];
}

int main(void)
{
  double h=0.71, Om=0.27, OL=0.73, Ob=0.046, XH=0.76;
  double H0_s = h*100.0*1e5/3.0857e24;
  double rho_crit0 = 1.8788e-29*h*h;
  double nH0 = Ob*rho_crit0*XH/mh;
  load_recfast("/Users/tabel/Projects/cicass/vbc_transfer/recfast/xeTrecfast.out");

  code_units u;
  u.comoving_coordinates = 0;
  u.density_units = mh;
  u.length_units  = 3.0857e21;
  u.time_units    = 3.1557e13;
  u.a_units       = 1.0;
  double z = 1000.0;
  u.a_value = 1.0/(1.0+z);
  set_velocity_units(&u);

  chemistry_data *cd = malloc(sizeof(chemistry_data));
  set_default_chemistry_parameters(cd);
  grackle_data->use_grackle = 1;
  grackle_data->with_radiative_cooling = 1;
  grackle_data->primordial_chemistry = 1;     /* H, He, e : recombination only */
  grackle_data->metal_cooling = 0;
  grackle_data->UVbackground = 0;
  grackle_data->cmb_recombination = atoi(getenv("CMBREC")?getenv("CMBREC"):"1");
  grackle_data->cosmology_hubble_constant_now = 71.0;
  grackle_data->cosmology_omega_matter_now    = 0.27;
  grackle_data->cosmology_omega_lambda_now    = 0.73;
  grackle_data->grackle_data_file = "../../input/CloudyData_noUVB.h5";
  if (initialize_chemistry_data(&u) == 0) { fprintf(stderr,"init failed\n"); return 1; }

  grackle_field_data f; int dim[3]={1,1,1}, st[3]={0,0,0}, en[3]={0,0,0};
  f.grid_rank=3; f.grid_dimension=dim; f.grid_start=st; f.grid_end=en; f.grid_dx=0.0;
  gr_float dens[1],ie[1],vx[1]={0},vy[1]={0},vz[1]={0};
  gr_float HI[1],HII[1],HeI[1],HeII[1],HeIII[1],e[1];
  f.density=dens; f.internal_energy=ie; f.x_velocity=vx; f.y_velocity=vy; f.z_velocity=vz;
  f.HI_density=HI; f.HII_density=HII; f.HeI_density=HeI; f.HeII_density=HeII; f.HeIII_density=HeIII;
  f.e_density=e;
  f.HM_density=NULL; f.H2I_density=NULL; f.H2II_density=NULL;
  f.DI_density=NULL; f.DII_density=NULL; f.HDI_density=NULL; f.metal_density=NULL;
  f.volumetric_heating_rate=NULL; f.specific_heating_rate=NULL;
  f.RT_HI_ionization_rate=NULL; f.RT_HeI_ionization_rate=NULL; f.RT_HeII_ionization_rate=NULL;
  f.RT_H2_dissociation_rate=NULL; f.RT_heating_rate=NULL;

  /* RECFAST ICs at z=1000. */
  double T = interp_rf(1000.0,Tt), xe = interp_rf(1000.0,xt), tiny=1e-20;
  double nH = nH0*pow(1.0+z,3.0), rho_tot = nH*mh/XH;
  dens[0] = rho_tot/u.density_units;
  HI[0]=(1.0-xe)*XH*dens[0]; HII[0]=xe*XH*dens[0]; e[0]=xe*XH*dens[0];
  HeI[0]=(1.0-XH)*dens[0]; HeII[0]=tiny*dens[0]; HeIII[0]=tiny*dens[0];
  double Tunits=get_temperature_units(&u), mu=1.22;
  ie[0] = T/Tunits/(5.0/3.0-1.0)/mu;

  printf("# CMBREC=%d\n# z   x_e(grackle)  x_e(recfast)  ratio   T_g(grk)  T(rf)\n",
         grackle_data->cmb_recombination);
  double dlna=0.01;
  while (z > 100.0) {
    double znew=(1.0+z)*exp(-dlna)-1.0; if(znew<100.0) znew=100.0;
    double fac=(1.0+znew)/(1.0+z), f3=fac*fac*fac;
    dens[0]*=f3; HI[0]*=f3; HII[0]*=f3; HeI[0]*=f3; HeII[0]*=f3; HeIII[0]*=f3; e[0]*=f3;
    ie[0]*=fac*fac;
    u.a_value=1.0/(1.0+znew);
    double dt=(cosmic_time(znew,H0_s,Om)-cosmic_time(z,H0_s,Om))/u.time_units;
    if (solve_chemistry(&u,&f,dt)==0){fprintf(stderr,"solve failed z=%g\n",z);return 1;}
    double nHn=HI[0]+HII[0];
    double xeg=e[0]/nHn;                       /* n_e/n_H (e carries m_H) */
    gr_float Tg[1]; calculate_temperature(&u,&f,Tg);
    double xrf=interp_rf(znew,xt)/XH;          /* table x_e is per H nucleus already? */
    /* recfast x_e is n_e/n_H; our xeg = e/(HI+HII) = n_e/n_H too. compare directly. */
    xrf=interp_rf(znew,xt);
    printf("%7.2f  %.4e  %.4e  %.4f  %8.2f  %8.2f\n",
           znew, xeg, xrf, xeg/xrf, (double)Tg[0], interp_rf(znew,Tt));
    z=znew;
  }
  return 0;
}
