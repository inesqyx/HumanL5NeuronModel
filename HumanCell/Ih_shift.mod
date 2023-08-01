  :Comment :
:Reference : :		Kole,Hallermann,and Stuart, J. Neurosci. 2006

NEURON	{
	SUFFIX Ih_shift
	NONSPECIFIC_CURRENT ihcn
	RANGE gIhbar, gIh, ihcn 
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gIhbar = 5.135e-05 (S/cm2) 
	ehcn =  -49.846 (mV)
	vh= -90.873
	k= 8.0488
	a= 23.453
	b= 0.21851
	c= 1.3091e-09
	d= 0.082671
	e= 1.4975e-09
	shift_minf = 0
	shift_mtau = 0
}


ASSIGNED	{
	v	(mV)
	ihcn	(mA/cm2)
	gIh	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gIh = gIhbar*m
	ihcn = gIh*(v-ehcn)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates(){
	UNITSOFF
		mInf=1/(1+exp(((v+shift_minf)-vh)/k))
		mTau=1/(exp(-a-b*(v+shift_mtau))+exp(-c+d*(v+shift_mtau)))+e
	UNITSON
}
