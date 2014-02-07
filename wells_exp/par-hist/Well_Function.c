
double Well_Func_u ( double t, void *params)
{
	struct u_params *p = (struct u_params *) params;

	double T1 = p->T1;
	double T2 = p->T2;
	double S1 = p->S1;
	double S2 = p->S2;
	double r = p->r;

	return pow(r,2) * S1 * exp(S2 / t)/( 4 * T1 * exp(T2 / t) * t);
}


