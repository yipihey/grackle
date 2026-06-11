/* Verify the neutral_helium T>1e4 K one-time warning fires. */
#include <stdlib.h>
#include <stdio.h>
#include <grackle.h>
#define mh 1.67262171e-24
int main(void){
  code_units u; u.comoving_coordinates=0; u.density_units=mh;
  u.length_units=3.0857e21; u.time_units=3.1557e13; u.a_units=1.0; u.a_value=1.0;
  set_velocity_units(&u);
  chemistry_data *cd=malloc(sizeof(chemistry_data));
  set_default_chemistry_parameters(cd);
  grackle_data->use_grackle=1; grackle_data->with_radiative_cooling=1;
  grackle_data->primordial_chemistry=2; grackle_data->UVbackground=0;
  grackle_data->neutral_helium=1;
  grackle_data->grackle_data_file="../../input/CloudyData_noUVB.h5";
  if(initialize_chemistry_data(&u)==0){fprintf(stderr,"init failed\n");return 1;}
  grackle_field_data f; int dim[3]={1,1,1},st[3]={0,0,0},en[3]={0,0,0};
  f.grid_rank=3; f.grid_dimension=dim; f.grid_start=st; f.grid_end=en; f.grid_dx=0;
  gr_float dens[1]={1.0},ie[1],vx[1]={0},vy[1]={0},vz[1]={0};
  gr_float HI[1],HII[1],HeI[1],HeII[1],HeIII[1],e[1],HM[1],H2I[1],H2II[1];
  f.density=dens; f.internal_energy=ie; f.x_velocity=vx; f.y_velocity=vy; f.z_velocity=vz;
  f.HI_density=HI; f.HII_density=HII; f.HeI_density=HeI; f.HeII_density=HeII;
  f.HeIII_density=HeIII; f.e_density=e; f.HM_density=HM; f.H2I_density=H2I;
  f.H2II_density=H2II; f.metal_density=NULL;
  f.volumetric_heating_rate=NULL; f.specific_heating_rate=NULL;
  f.RT_HI_ionization_rate=NULL; f.RT_HeI_ionization_rate=NULL;
  f.RT_HeII_ionization_rate=NULL; f.RT_H2_dissociation_rate=NULL; f.RT_heating_rate=NULL;
  double XH=0.76,t=1e-20;
  HI[0]=0.5*XH; HII[0]=0.5*XH; e[0]=0.5*XH; HeI[0]=(1-XH); HeII[0]=t; HeIII[0]=t;
  HM[0]=t; H2I[0]=t; H2II[0]=t;
  /* T = 2e4 K (above the 1e4 validity ceiling) */
  double Tu=get_temperature_units(&u);
  ie[0]= 2.0e4 / Tu / (5.0/3.0-1.0) / 1.0;
  fprintf(stderr,"-- calling solve_chemistry at T=2e4 K with neutral_helium=1 --\n");
  solve_chemistry(&u,&f,1e-3);
  fprintf(stderr,"-- done --\n");
  return 0;
}
