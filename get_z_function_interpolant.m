function [Z_interp] = get_z_function_interpolant()

z_function_table = 'zFunction.nc';

if is_octave()
    table_zeta_re = ncread(z_function_table,'arg_re');
    table_Z_re    = ncread(z_function_table,'Z_re');
    table_Z_im    = ncread(z_function_table,'Z_im');
    table_Zp_re   = ncread(z_function_table,'Zp_re');
    table_Zp_im   = ncread(z_function_table,'Zp_im');
    
    pp_Z_re  = spline(table_zeta_re,table_Z_re);
    pp_Z_im  = spline(table_zeta_re,table_Z_im);
    pp_Zp_re = spline(table_zeta_re,table_Zp_re);
    pp_Zp_im = spline(table_zeta_re,table_Zp_im);
    
    Z_interp.Z_re  = pp_Z_re;
    Z_interp.Z_im  = pp_Z_im;
    Z_interp.Zp_re = pp_Zp_re;
    Z_interp.Zp_im = pp_Zp_im;
    
    Z_interp.table_arg_re = table_zeta_re;
    Z_interp.table_Z_re  = table_Z_re;
    Z_interp.table_Z_im  = table_Z_im;
    Z_interp.table_Zp_re = table_Zp_re;
    Z_interp.table_Zp_im = table_Zp_im;
    
else
    ncid = netcdf.open(z_function_table);
    
    arg_re_id = netcdf.inqVarID(ncid,'arg_re');
    Z_re_id = netcdf.inqVarID(ncid,'Z_re');
    Z_im_id = netcdf.inqVarID(ncid,'Z_im');
    Zp_re_id = netcdf.inqVarID(ncid,'Zp_re');
    Zp_im_id = netcdf.inqVarID(ncid,'Zp_im');
    
    table_zeta_re = netcdf.getVar(ncid,arg_re_id);
    table_Z_re = netcdf.getVar(ncid,Z_re_id);
    table_Z_im = netcdf.getVar(ncid,Z_im_id);
    table_Zp_re = netcdf.getVar(ncid,Zp_re_id);
    table_Zp_im = netcdf.getVar(ncid,Zp_im_id);
    
    int_Z_re = griddedInterpolant(table_zeta_re,table_Z_re,'spline');
    int_Z_im = griddedInterpolant(table_zeta_re,table_Z_im,'spline');
    int_Zp_re = griddedInterpolant(table_zeta_re,table_Zp_re,'spline');
    int_Zp_im = griddedInterpolant(table_zeta_re,table_Zp_im,'spline');
    
    Z_interp.Z_re = int_Z_re;
    Z_interp.Z_im = int_Z_im;
    Z_interp.Zp_re = int_Zp_re;
    Z_interp.Zp_im = int_Zp_im;
    
end

min_zeta_re = min(table_zeta_re);
max_zeta_re = max(table_zeta_re);

min_zeta_re = min(table_zeta_re);
max_zeta_re = max(table_zeta_re);

Z_interp.min_zeta_re = min_zeta_re;
Z_interp.max_zeta_re = max_zeta_re;

end