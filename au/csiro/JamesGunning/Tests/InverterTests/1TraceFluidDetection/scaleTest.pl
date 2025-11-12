#!/usr/bin/perl
# Script to examine dependence of solver on scaling
for ($p=-2;$p<=2;$p+=0.2)
{
	$scale=exp($p*log(10.0));
	print $p," ",$scale, " ";
	system("java au.csiro.JamesGunning.seismicinverter.SeismicInverter ModelDescription.xml RockandFluidProperties.xml -v 4 -m prior_traces.su -N 1 -PD -IF --Dval $scale -D1 --secant -II| awk '/Minimiser its =/ &&p==0{print \$4;p=1}'");
}
