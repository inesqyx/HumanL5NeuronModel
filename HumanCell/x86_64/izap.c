/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__Izap
#define _nrn_initial _nrn_initial__Izap
#define nrn_cur _nrn_cur__Izap
#define _nrn_current _nrn_current__Izap
#define nrn_jacob _nrn_jacob__Izap
#define nrn_state _nrn_state__Izap
#define _net_receive _net_receive__Izap 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define del _p[0]
#define dur _p[1]
#define f0 _p[2]
#define f1 _p[3]
#define amp _p[4]
#define f _p[5]
#define i _p[6]
#define on _p[7]
#define v _p[8]
#define _g _p[9]
#define _tsav _p[10]
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "del", "ms",
 "dur", "ms",
 "f0", "1/s",
 "f1", "1/s",
 "amp", "nA",
 "f", "1/s",
 "i", "nA",
 0,0
};
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Izap",
 "del",
 "dur",
 "f0",
 "f1",
 "amp",
 0,
 "f",
 "i",
 0,
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 11, _prop);
 	/*initialize range parameters*/
 	del = 0;
 	dur = 0;
 	f0 = 0;
 	f1 = 0;
 	amp = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 11;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 
#define _tqitem &(_ppvar[2]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _izap_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 11, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Izap /Users/scottrich/OneDrive - UHN/3) NEURON stuff/HumanCellOrganized_v2/x86_64/izap.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double PI = 3.14159;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{  double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _thread = (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 1.0 ) {
     if ( on  == 0.0 ) {
       on = 1.0 ;
       net_send ( _tqitem, _args, _pnt, t +  dur , 1.0 ) ;
       }
     else {
       on = 0.0 ;
       }
     }
   } }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
 {
   f = 0.0 ;
   i = 0.0 ;
   on = 0.0 ;
   if ( del < 0.0 ) {
     del = 0.0 ;
     }
   if ( dur < 0.0 ) {
     dur = 0.0 ;
     }
   if ( f0 <= 0.0 ) {
     f0 = 0.0 ;
     }
   if ( f1 <= 0.0 ) {
     f1 = 0.0 ;
     }
   if ( dur > 0.0 ) {
     net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  del , 1.0 ) ;
     }
   }

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   if ( on  == 0.0 ) {
     f = 0.0 ;
     i = 0.0 ;
     }
   else {
     f = f0 + ( f1 - f0 ) * ( t - del ) / dur ;
     i = amp * sin ( 2.0 * PI * ( t - del ) * ( f0 + ( f1 - f0 ) * ( t - del ) / ( 2.0 * dur ) ) * ( 0.001 ) ) ;
     }
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) += _rhs;
  }else
#endif
  {
	NODERHS(_nd) += _rhs;
  }
  if (_nt->_nrn_fast_imem) { _nt->_nrn_fast_imem->_nrn_sav_rhs[_ni[_iml]] += _rhs; }
#if EXTRACELLULAR
 if (_nd->_extnode) {
   *_nd->_extnode->_rhs[0] += _rhs;
 }
#endif
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) -= _g;
  }else
#endif
  {
	NODED(_nd) -= _g;
  }
  if (_nt->_nrn_fast_imem) { _nt->_nrn_fast_imem->_nrn_sav_d[_ni[_iml]] -= _g; }
#if EXTRACELLULAR
 if (_nd->_extnode) {
   *_nd->_extnode->_d[0] += _g;
 }
#endif
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/scottrich/OneDrive - UHN/3) NEURON stuff/HumanCellOrganized_v2/izap.mod";
static const char* nmodl_file_text = 
  ": $Id: izap.mod,v 1.4 2015/05/28 02:42:00 ted Exp ted $\n"
  "\n"
  "COMMENT\n"
  "izap.mod\n"
  "\n"
  "Delivers an oscillating current that starts at t = del >= 0.\n"
  "The frequency of the oscillation increases linearly with time\n"
  "from f0 at t == del to f1 at t == del + dur, \n"
  "where both del and dur are > 0.\n"
  "\n"
  "Uses event delivery system to ensure compatibility with adaptive integration.\n"
  "\n"
  "Original implementation 12/4/2008 NTC\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "  POINT_PROCESS Izap\n"
  "  RANGE del, dur, f0, f1, amp, f, i\n"
  "  ELECTRODE_CURRENT i\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "  (nA) = (nanoamp)\n"
  "  PI = (pi) (1)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "  del (ms)\n"
  "  dur (ms)\n"
  "  f0 (1/s)  : frequency is in Hz\n"
  "  f1 (1/s)\n"
  "  amp (nA)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "  f (1/s)\n"
  "  i (nA)\n"
  "  on (1)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "  f = 0\n"
  "  i = 0\n"
  "  on = 0\n"
  "\n"
  "  if (del<0) { del=0 }\n"
  "  if (dur<0) { dur=0 }\n"
  "  if (f0<=0) { f0=0 (1/s) }\n"
  "  if (f1<=0) { f1=0 (1/s) }\n"
  "\n"
  "  : do nothing if dur == 0\n"
  "  if (dur>0) {\n"
  "    net_send(del, 1)  : to turn it on and start frequency ramp\n"
  "  }\n"
  "}\n"
  "\n"
  "COMMENT\n"
  "The angular velocity in radians/sec is w = 2*PI*f, \n"
  "where f is the instantaneous frequency in Hz.\n"
  "\n"
  "Assume for the moment that the frequency ramp starts at t = 0.\n"
  "f = f0 + (f1 - f0)*t/dur\n"
  "\n"
  "Then the angular displacement is\n"
  "theta = 2*PI * ( f0*t + (f1 - f0)*(t^2)/(2*dur) ) \n"
  "      = 2*PI * t * (f0 + (f1 - f0)*t/(2*dur))\n"
  "But the ramp starts at t = del, so just substitute t-del for every occurrence of t\n"
  "in the formula for theta.\n"
  "ENDCOMMENT\n"
  "\n"
  "BREAKPOINT {\n"
  "  if (on==0) {\n"
  "    f = 0\n"
  "    i = 0\n"
  "  } else {\n"
  "    f = f0 + (f1 - f0)*(t-del)/dur\n"
  "    i = amp * sin( 2*PI * (t-del) * (f0 + (f1 - f0)*(t-del)/(2*dur)) * (0.001) )\n"
  "  }\n"
  "}\n"
  "\n"
  "NET_RECEIVE (w) {\n"
  "  : respond only to self-events with flag > 0\n"
  "  if (flag == 1) {\n"
  "    if (on==0) {\n"
  "      on = 1  : turn it on\n"
  "      net_send(dur, 1)  : to stop frequency ramp, freezing frequency at f1\n"
  "    } else {\n"
  "      on = 0  : turn it off\n"
  "    }\n"
  "  }\n"
  "}\n"
  ;
#endif
