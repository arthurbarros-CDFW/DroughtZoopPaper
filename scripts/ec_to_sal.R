#EC to salinity ppt
ec_2_sal = function(temp, cond){
  if (any(temp > 35, na.rm = T)) {
    warning('Temperature is high, ensure that units are in degrees C', call. = F, immediate. = T)
  }
  if(length(temp == 1)) {
    temp = rep(temp, length.out = length(cond))
  } else if ( length(temp == length(cond))) {
    temp = temp
  } else {
    warning('Temperature must be of length 1 or match number of conductivity readings')
  }
  ref_cond = 42914
  cond_rat = cond/ref_cond
  rt = 0.6766097 + (0.0200564*temp) + (0.0001104259*(temp^2)) + ((-6.9698*10^-7)*(temp^3)) + ((1.0031*10^-9)*temp^4)
  Rt = cond_rat/rt
  dS = ((temp-15)/(1+0.0162*(temp-15)))*(0.0005+(-0.0056)*(Rt^0.5)+(-0.0066)*Rt+(-0.0375)*(Rt^1.5)+(0.0636)*(Rt^2)+(-0.0144)*(Rt^2.5))
  sal = c()
  sal[is.na(cond)] = NA
  sal[cond >3000 & !is.na(cond)] = 0.008 + (-0.1692)*(Rt[cond>3000 & !is.na(cond)]^0.5) + 25.3851*Rt[cond > 3000 & !is.na(cond)]+14.0941*(Rt[cond > 3000 & !is.na(cond)]^1.5)+(-7.0261)*(Rt[cond > 3000 & !is.na(cond)]^2)+2.7081*(Rt[cond > 3000 & !is.na(cond)]^2.5)+dS[cond > 3000 & !is.na(cond)]
  sal[cond <=3000 & !is.na(cond)] = (0.008+ (-0.1692)*(Rt[cond <= 3000 & !is.na(cond)]^0.5)+25.3851*Rt[cond <= 3000 & !is.na(cond)]+14.0941*(Rt[cond <= 3000 & !is.na(cond)]^1.5)+(-7.0261)*(Rt[cond <= 3000 & !is.na(cond)]^2)+2.7081*(Rt[cond <= 3000 & !is.na(cond)]^2.5)+dS[cond <= 3000 & !is.na(cond)]) -
    (0.008/(1+(1.5*(400*Rt[cond <= 3000 & !is.na(cond)]))+((400*Rt[cond <= 3000 & !is.na(cond)])^2))-(0.0005*(temp[cond <= 3000 & !is.na(cond)]-15)/(1+0.0162*(temp[cond <= 3000 & !is.na(cond)]-15)))/(1+((100*Rt[cond <= 3000 & !is.na(cond)])^0.5)+((100*Rt[cond <= 3000 & !is.na(cond)])^1.5)))
  return(sal)
}
