/* One-zone validation of the CMB photo-destruction of H- and H2+.
   A post-recombination primordial parcel (z fixed, no UV background) is evolved
   with cmb_dissociation off (argv1=0) then on (argv1=1). With the CMB rates on,
   the H- and H2+ intermediaries are photo-destroyed and H2 formation is
   suppressed. Usage: ./cmb_test <0|1> <redshift>
   Field allocation mirrors c_example.c (all fields allocated). */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grackle.h>

#define mh 1.67262171e-24

int main(int argc, char *argv[])
{
  int cmb = (argc > 1) ? atoi(argv[1]) : 0;
  double zred = (argc > 2) ? atof(argv[2]) : 200.0;
  grackle_verbose = 0;

  code_units u;
  u.comoving_coordinates = 0;
  u.density_units = 0.26 * mh;          /* ~ mean baryon density at z=200 */
  u.length_units  = 1.0;
  u.time_units    = 1.0e12;
  u.a_units       = 1.0;
  u.a_value       = 1.0 / (1.0 + zred); /* carries z -> T_CMB = 2.73*(1+z) */
  set_velocity_units(&u);

  chemistry_data *cd = malloc(sizeof(chemistry_data));
  set_default_chemistry_parameters(cd);
  grackle_data->use_grackle = 1;
  grackle_data->with_radiative_cooling = 1;
  grackle_data->primordial_chemistry = 2;
  grackle_data->metal_cooling = 0;
  grackle_data->dust_chemistry = 0;
  grackle_data->UVbackground = 0;
  grackle_data->cmb_dissociation = cmb;       /* <-- the new flag */
  grackle_data->grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5";
  if (initialize_chemistry_data(&u) == 0) { fprintf(stderr,"init failed\n"); return 1; }

  double tiny = 1e-20;
  grackle_field_data f;
  int dim[3], st[3], en[3];
  f.grid_rank = 3; f.grid_dimension = dim; f.grid_start = st; f.grid_end = en;
  f.grid_dx = 0.0;
  for (int i=0;i<3;i++){ dim[i]=1; st[i]=0; en[i]=0; }

  int N = 1;
  f.density         = malloc(N*sizeof(gr_float));
  f.internal_energy = malloc(N*sizeof(gr_float));
  f.x_velocity      = malloc(N*sizeof(gr_float));
  f.y_velocity      = malloc(N*sizeof(gr_float));
  f.z_velocity      = malloc(N*sizeof(gr_float));
  f.HI_density      = malloc(N*sizeof(gr_float));
  f.HII_density     = malloc(N*sizeof(gr_float));
  f.HeI_density     = malloc(N*sizeof(gr_float));
  f.HeII_density    = malloc(N*sizeof(gr_float));
  f.HeIII_density   = malloc(N*sizeof(gr_float));
  f.e_density       = malloc(N*sizeof(gr_float));
  f.HM_density      = malloc(N*sizeof(gr_float));
  f.H2I_density     = malloc(N*sizeof(gr_float));
  f.H2II_density    = malloc(N*sizeof(gr_float));
  f.DI_density      = malloc(N*sizeof(gr_float));
  f.DII_density     = malloc(N*sizeof(gr_float));
  f.HDI_density     = malloc(N*sizeof(gr_float));
  f.metal_density   = malloc(N*sizeof(gr_float));
  f.volumetric_heating_rate = malloc(N*sizeof(gr_float));
  f.specific_heating_rate   = malloc(N*sizeof(gr_float));
  f.RT_HI_ionization_rate   = malloc(N*sizeof(gr_float));
  f.RT_HeI_ionization_rate  = malloc(N*sizeof(gr_float));
  f.RT_HeII_ionization_rate = malloc(N*sizeof(gr_float));
  f.RT_H2_dissociation_rate = malloc(N*sizeof(gr_float));
  f.RT_heating_rate         = malloc(N*sizeof(gr_float));

  double Tunits = get_temperature_units(&u);
  double xH = grackle_data->HydrogenFractionByMass;
  double xe0 = 2.0e-4;                          /* residual electron fraction */
  f.density[0] = 1.0;
  f.HII_density[0] = xe0*xH;     f.e_density[0] = xe0*xH;
  f.HI_density[0]  = xH - xe0*xH;
  f.HeI_density[0] = 1.0 - xH;   f.HeII_density[0] = tiny; f.HeIII_density[0] = tiny;
  f.HM_density[0]  = tiny; f.H2I_density[0] = tiny; f.H2II_density[0] = tiny;
  f.DI_density[0]  = tiny; f.DII_density[0] = tiny; f.HDI_density[0] = tiny;
  f.metal_density[0] = tiny;
  f.x_velocity[0]=0; f.y_velocity[0]=0; f.z_velocity[0]=0;
  f.internal_energy[0] = (2.73*(1.0+zred)) / Tunits / (5.0/3.0 - 1.0) / 1.22; /* T = T_CMB */
  f.volumetric_heating_rate[0]=0; f.specific_heating_rate[0]=0;
  f.RT_HI_ionization_rate[0]=0; f.RT_HeI_ionization_rate[0]=0;
  f.RT_HeII_ionization_rate[0]=0; f.RT_H2_dissociation_rate[0]=0; f.RT_heating_rate[0]=0;

  double dt = 3.15e7 * 5e4 / u.time_units;       /* 5e4 yr/step */
  for (int n = 0; n < 1000; n++)
    if (solve_chemistry(&u, &f, dt) == 0) { fprintf(stderr,"solve failed\n"); return 1; }

  printf("cmb=%d z=%.0f  H2I/rho=%.4e  HM/rho=%.4e  e/rho=%.4e\n",
         cmb, zred, (double)f.H2I_density[0]/f.density[0],
         (double)f.HM_density[0]/f.density[0], (double)f.e_density[0]/f.density[0]);
  return 0;
}
