/* Homogeneous primordial-chemistry evolution through cosmic expansion, z=1000 -> 10,
   using the Grackle network (the same chemistry Enzo uses) with our CMB photo-rates.
   Reproduces the Galli & Palla (1998) species-fraction-vs-redshift figure.

   Grackle solves the chemistry + Compton + radiative cooling each step; the caller
   advances the expansion: physical densities dilute (x (1+z)^3) and the gas cools
   adiabatically (T x (1+z)^2) between steps. Output: a table of number fractions
   n(X)/n_H vs z.  Usage: ./cosmo_evolve > cosmo_species.dat */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grackle.h>

#define mh 1.67262171e-24
#define kb 1.3806504e-16

/* EdS cosmic time [s] (matter-dominated, accurate for z>10). */
static double cosmic_time(double z, double H0_s, double Om) {
  return (2.0/3.0)/H0_s/sqrt(Om) * pow(1.0+z, -1.5);
}

int main(void)
{
  /* cosmology */
  double h = 0.71, Om = 0.27, Ob = 0.046, XH = 0.76, DtoH = 2.6e-5;
  double H0_s = h*100.0*1e5/3.0857e24;                 /* H0 [1/s] */
  double rho_crit0 = 1.8788e-29*h*h;                   /* g/cm^3 */
  double nH0 = Ob*rho_crit0*XH/mh;                     /* mean n_H today [cm^-3] */

  /* units (physical; a_value carries z) */
  code_units u;
  u.comoving_coordinates = 0;
  u.density_units = mh;            /* density field ~ number density of H-mass units */
  u.length_units  = 3.0857e21;     /* 1 kpc  (length/time -> physical velocity unit) */
  u.time_units    = 3.1557e13;     /* 1 Myr  -> velocity_units ~ 978 km/s, T_units ~ 1e8 K */
  u.a_units       = 1.0;
  double z = 1000.0;
  u.a_value = 1.0/(1.0+z);
  set_velocity_units(&u);

  chemistry_data *cd = malloc(sizeof(chemistry_data));
  set_default_chemistry_parameters(cd);
  grackle_data->use_grackle = 1;
  grackle_data->with_radiative_cooling = 1;
  grackle_data->primordial_chemistry = atoi(getenv("PC")?getenv("PC"):"3"); /* PC=1 -> no molecules (cooling diagnostic) */
  grackle_data->metal_cooling = 0;
  grackle_data->dust_chemistry = 0;
  grackle_data->UVbackground = 0;
  grackle_data->cmb_dissociation = atoi(getenv("CMBDISS")?getenv("CMBDISS"):"1"); /* CMB H-/H2+ photo-destruction */
  grackle_data->equilibrium_h2_intermediates = atoi(getenv("EQUIL")?getenv("EQUIL"):"0"); /* v2026 */
  grackle_data->cmb_recombination = atoi(getenv("CMBREC")?getenv("CMBREC"):"1"); /* Peebles C-factor */
  grackle_data->neutral_helium = atoi(getenv("HE")?getenv("HE"):"0"); /* HE=1: He all neutral, n_e=n_HII */
  grackle_data->equilibrium_deuterium = atoi(getenv("DEUT")?getenv("DEUT"):"0"); /* DEUT=1: advect only HD */
  grackle_data->cosmology_hubble_constant_now = 71.0;   /* km/s/Mpc (for H(z) in the C-factor) */
  grackle_data->cosmology_omega_matter_now    = 0.27;
  grackle_data->cosmology_omega_lambda_now    = 0.73;
  grackle_data->grackle_data_file = "../../input/CloudyData_noUVB.h5";
  if (initialize_chemistry_data(&u) == 0) { fprintf(stderr,"init failed\n"); return 1; }

  /* one zone */
  grackle_field_data f; int dim[3]={1,1,1}, st[3]={0,0,0}, en[3]={0,0,0};
  f.grid_rank=3; f.grid_dimension=dim; f.grid_start=st; f.grid_end=en; f.grid_dx=0.0;
  gr_float dens[1]; f.density=dens;
  gr_float ie[1];   f.internal_energy=ie;
  gr_float vx[1]={0},vy[1]={0},vz[1]={0}; f.x_velocity=vx; f.y_velocity=vy; f.z_velocity=vz;
  gr_float HI[1],HII[1],HeI[1],HeII[1],HeIII[1],e[1],HM[1],H2I[1],H2II[1],DI[1],DII[1],HDI[1],Z[1];
  f.HI_density=HI; f.HII_density=HII; f.HeI_density=HeI; f.HeII_density=HeII; f.HeIII_density=HeIII;
  f.e_density=e; f.HM_density=HM; f.H2I_density=H2I; f.H2II_density=H2II;
  f.DI_density=DI; f.DII_density=DII; f.HDI_density=HDI; f.metal_density=Z;
  f.volumetric_heating_rate=NULL; f.specific_heating_rate=NULL;
  f.RT_HI_ionization_rate=NULL; f.RT_HeI_ionization_rate=NULL; f.RT_HeII_ionization_rate=NULL;
  f.RT_H2_dissociation_rate=NULL; f.RT_heating_rate=NULL;

  /* initial state at z=1000 from RECFAST: x_e=0.047, T=2728 K.
     comoving_coordinates=0 -> the density field is PHYSICAL (grackle does not scale
     by a^3; a_value only sets the redshift / CMB temperature). We dilute the physical
     density manually as the universe expands (x (1+z)^3 -> x fac^3 per step). */
  double T = 2728.0, xe = 0.047, tiny = 1e-20;
  double nH = nH0*pow(1.0+z, 3.0);                      /* physical n_H at z=1000 [cm^-3] */
  double rho_tot = nH*mh/XH;                            /* physical total mass incl He */
  dens[0] = rho_tot/u.density_units;
  HI[0]   = (1.0-xe)*XH*dens[0];  HII[0] = xe*XH*dens[0];  e[0] = xe*XH*dens[0];
  HeI[0]  = (1.0-XH)*dens[0];     HeII[0]=tiny*dens[0]; HeIII[0]=tiny*dens[0];
  /* seed x_H2 = n(H2)/n_H = 1e-15  (H2I is the H2 mass density = 2*n(H2)). */
  HM[0]=tiny*dens[0]; H2I[0]=2.0e-15*XH*dens[0]; H2II[0]=tiny*dens[0];
  DI[0]=(1.0-xe)*DtoH*XH*dens[0]; DII[0]=xe*DtoH*XH*dens[0]; HDI[0]=tiny*dens[0]; Z[0]=tiny*dens[0];
  double Tunits = get_temperature_units(&u), mu = 1.22;
  ie[0] = T / Tunits / (5.0/3.0-1.0) / mu;

  printf("# z T_gas T_cmb x_HII x_H2 x_HM x_H2p x_e x_HD x_DII\n");
  double dlna = 0.004;                                  /* expansion step */
  while (z > 10.0) {
    double znew = (1.0+z)*exp(-dlna) - 1.0;
    if (znew < 10.0) znew = 10.0;
    double fac = (1.0+znew)/(1.0+z);                    /* < 1 */
    /* expansion: dilute physical densities (x fac^3), adiabatic cool (T x fac^2) */
    double f3 = fac*fac*fac;
    dens[0]*=f3; HI[0]*=f3; HII[0]*=f3; HeI[0]*=f3; HeII[0]*=f3; HeIII[0]*=f3; e[0]*=f3;
    HM[0]*=f3; H2I[0]*=f3; H2II[0]*=f3; DI[0]*=f3; DII[0]*=f3; HDI[0]*=f3; Z[0]*=f3;
    ie[0] *= fac*fac;
    u.a_value = 1.0/(1.0+znew);
    double dt = (cosmic_time(znew,H0_s,Om) - cosmic_time(z,H0_s,Om))/u.time_units;
    if (solve_chemistry(&u, &f, dt) == 0) { fprintf(stderr,"solve failed z=%g\n",z); return 1; }

    /* number fractions relative to n_H (grackle: H2/He/HD carry mass; e counts m_H) */
    double nHn = HI[0]+HII[0]+HM[0]+H2I[0]+H2II[0];     /* H nuclei (field units) */
    gr_float Tg[1]; calculate_temperature(&u,&f,Tg);
    printf("%.4f %.3e %.3e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
           znew, (double)Tg[0], 2.73*(1.0+znew),
           HII[0]/nHn, (H2I[0]/2.0)/nHn, HM[0]/nHn, (H2II[0]/2.0)/nHn,
           e[0]/nHn, (HDI[0]/3.0)/nHn, DII[0]/nHn);
    z = znew;
  }
  return 0;
}
