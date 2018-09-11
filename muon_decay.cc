#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "LHAPDF/LHAPDF.h"

#include <boost/progress.hpp>

#include "TRandom.h"
#include "TRandom3.h"
#include "TLorentzVector.h" 
#include "TVector3.h"

using namespace ThePEG;
using namespace LHAPDF;

double sqr(double a){ return a*a; }

// solve v = (n+2) * u^(n+1) - (n+1) * u^(n+2) for u
double bisect(double v, double n, 
	      double target = -16., double maxLevel = 80.) {

  if ( v != 0.0 && v != 1.0 ) {

    double level = 0;
    double left = 0;
    double right = 1;

    double checkV = -1.;
    double u = -1;

    while ( level < maxLevel ) {

      u = (left+right)*pow(0.5,level+1.);
      checkV = 
	pow(u,n+1.)*(n+2.-(n+1.)*u);

      if ( log10(abs(1.-checkV/v)) <= target )
	break;

      left *= 2.;
      right *= 2.;

      if ( v <= checkV ) {
	right -= 1.;
	++level;
      }

      if ( v > checkV ) {
	left += 1.;
	++level;
      }

    }

    return u;

  }

  return v;

}

// solve the reshuffling equation
double bisectReshuffling(const vector<Lorentz5Momentum>& momenta,
			 const vector<Energy>& masses,
			 Energy w,
			 double target = -16., double maxLevel = 80.) {
  
  double level = 0;
  double left = 0;
  double right = 1;

  double check = -1.;
  double xi = -1;

  while ( level < maxLevel ) {

    xi = (left+right)*pow(0.5,level+1.);
    check = 0.;
    vector<Lorentz5Momentum>::const_iterator p = momenta.begin();
    vector<Energy>::const_iterator  m = masses.begin();
    for ( ; p != momenta.end(); ++p, ++m ) {
      check += sqrt(sqr(xi)*p->vect().mag2()+sqr(*m))/w;
    }

    if ( log10(abs(1.-check)) <= target )
      break;

    left *= 2.;
    right *= 2.;

    if ( check >= 1. ) {
      right -= 1.;
      ++level;
    }

    if ( check < 1. ) {
      left += 1.;
      ++level;
    }

  }

  return xi;  

}


// massless weights
static double weights[45] = {

  -1.,-1.,
  0.039788735772973833942,    // 2
  0.00012598255637968550463,  // 3
  1.3296564302788840628E-7,   // 4
  7.0167897579949011130E-11,  // 5
  2.2217170114046130768E-14,  // 6
  4.68973E-18, // 7
  7.07097E-22, // 8
  7.99597E-26, // 9
  7.03265E-30, // 10
  4.94831E-34, // 11
  2.84868E-38, // 12
  1.36663E-42, // 13
  5.54761E-47, // 14
  1.93026E-51, // 15
  5.82071E-56, // 16
  1.53584E-60, // 17
  3.57566E-65, // 18
  7.39972E-70, // 19
  1.37015E-74, // 20
  2.28332E-79, // 21
  3.44268E-84, // 22
  4.71884E-89, // 23
  5.90561E-94, // 24
  6.77495E-99, // 25
  7.15048E-104, // 26
  6.9663E-109, // 27
  6.28413E-114, // 28
  5.26385E-119, // 29
  4.10514E-124, // 30
  2.98806E-129, // 31
  2.03463E-134, // 32
  1.29884E-139, // 33
  7.78882E-145, // 34
  4.39601E-150, // 35
  2.33933E-155, // 36
  1.17571E-160, // 37
  5.58956E-166, // 38
  2.51752E-171, // 39
  1.07573E-176, // 40
  4.36677E-182, // 41
  1.68615E-187, // 42
  6.20075E-193, // 43
  2.17424E-198   // 44
};



// RAMBO generation
pair<double, vector<Lorentz5Momentum>> generateRamboKinematics(
             vector<Energy> m,
			       Energy Ecm) {

  vector<Lorentz5Momentum> P;
  for(int i=0; i<m.size(); i++){
    Lorentz5Momentum pdm;
    pdm.setMass(m[i]);
    P.push_back(pdm);
  }

  vector<double> r;
  r.resize(4*m.size());

  for ( vector<double>::iterator rnd = r.begin(); rnd != r.end(); ++rnd ) {
    *rnd = drand48();
  }

  Energy w = Ecm;
  size_t count = 0;
  Lorentz5Momentum Q;
  for ( vector<Lorentz5Momentum>::iterator k = P.begin();
	k != P.end(); ++k ) {
    Energy q = -w*log(r[count]*r[count+1]);
    double ct = 2.*r[count+2]-1.;
    double st = sqrt(1.-sqr(ct));
    double phi = 2.*Constants::pi*r[count+3];
    double cphi = cos(phi);
    double sphi = sqrt(1.-sqr(cphi));
    if ( phi > Constants::pi )
      sphi = -sphi;
    (*k).setMass(ZERO);
    (*k).setT(q);
    (*k).setX(q*cphi*st);
    (*k).setY(q*sphi*st);
    (*k).setZ(q*ct);
    count += 4;
    Q += *k;
  }

  Energy M = sqrt(Q.m2());
  double x = w/M;
  Boost beta = -(Q.vect() * (1./M));
  double gamma = Q.t()/M;
  double a = 1./(1.+gamma);

  for ( vector<Lorentz5Momentum>::iterator k = P.begin();
	k != P.end(); ++k ) {
    Energy q = (*k).t();
    Energy bq = beta*(*k).vect();
    (*k).setT(x*(gamma*q+bq));
    (*k).setVect(x*((*k).vect()+(q+a*bq)*beta));
  }

  size_t n = P.size();
  double weight = weights[n];

  double xi = bisectReshuffling(P,m,w);
  weight *= pow(xi,3.*(n-1.));

  Energy num = ZERO;
  Energy den = ZERO;

  vector<Energy>::const_iterator d = m.begin();
  for ( vector<Lorentz5Momentum>::iterator k = P.begin();
	k != P.end(); ++k, ++d ) {
    num += (*k).vect().mag2()/(*k).t();
    Energy q = (*k).t();
    (*k).setT(sqrt(sqr(*d)+xi*xi*sqr((*k).t())));
    (*k).setVect(xi*(*k).vect());
    weight *= q/(*k).t();
    den += (*k).vect().mag2()/(*k).t();
    (*k).setMass(*d);
    (*k).rescaleEnergy();
  }

  weight *= num/den;

  double scale = pow(Ecm/GeV, 2*n-4);
  double wei = weight*scale;
  //double wei = weight;

  pair<double, vector<Lorentz5Momentum>> res = make_pair(wei, P);
  //return weight*scale;
  return res;

}

double prod( Lorentz5Momentum p1, Lorentz5Momentum p2 ){
  double e1 = p1.t()/GeV;
  double x1 = p1.x()/GeV;
  double y1 = p1.y()/GeV;
  double z1 = p1.z()/GeV;
  double e2 = p2.t()/GeV;
  double x2 = p2.x()/GeV;
  double y2 = p2.y()/GeV;
  double z2 = p2.z()/GeV;
  return e1*e2 - x1*x2 - y1*y2 - z1*z2;
}

double ME2( Energy mmu, Lorentz5Momentum p2, Lorentz5Momentum q1, Lorentz5Momentum q2 ){
  double GF = 1.1663787 * pow(10., -5);
  double fac1 = 64.*GF*GF;
  double fac2 = mmu/GeV * q2.t()/GeV; 
  double fac3 = prod(p2, q1);  
  return fac1*fac2*fac3;
}

int main() {

  double pi = atan(1.)*4.;
  // size_t nf; cin >> nf;
  // size_t nw; cin >> nw;
  // double m; cin >> m;

  //#################################################################
  //  Muon decay
  //#################################################################
  int nev = 10000;

  Energy mmu = 0.1056583745*GeV;  
  vector<Energy> massvec;
  massvec.push_back(0.*GeV);
  massvec.push_back(0.*GeV);
  massvec.push_back(0.000*GeV);

  double GF = 1.1663787 * pow(10., -5);
  double fac1 = 64.*GF*GF;

  double mm = mmu/GeV;
  double dec_mu = GF*GF * pow(mm, 5) / (192.*pow(pi,3)); // 3.00918e-19 GeV

  double PS = 0;

  double dec = 0;

  for(int iev=0; iev < nev; iev++){  

    pair<double, vector<Lorentz5Momentum>> res = generateRamboKinematics(massvec, mmu);
    double wei = res.first;
    vector<Lorentz5Momentum> p = res.second;

    Lorentz5Momentum q1 = p[0];
    Lorentz5Momentum q2 = p[1];
    Lorentz5Momentum p2 = p[2];

    double me2_ = ME2(mmu, p2, q1, q2);
    //cout << wei <<"  "<<  me2_ << endl;

    dec += wei * me2_;

  }

  double dec_rambo = dec/nev;

  cout << "The muon decay rate validation:" << endl;
  cout << "[Rambo integration] = " << dec_rambo << endl;
  cout << "    [Exact formula] = " << dec_mu << endl;
  cout << "      [Rambo/Exact] = " << dec_rambo/dec_mu << endl;

  return 0;
  //#################################################################

  // int nev = 10000;

  // Energy Ecm = 9.*1000.*GeV;  
  // Energy mw = 80.4 * GeV;
  // int nf = 7;
  // int nw = 10;

  // vector<Energy> massvec;
  // for(int i=0; i<nf; i++) massvec.push_back(0.*GeV);
  // for(int i=0; i<nw; i++) massvec.push_back(mw);

  // double PS = 0;
  // for(int iev; iev < nev; iev++){  
  //   pair<double, vector<Lorentz5Momentum>> res = generateRamboKinematics(massvec, Ecm);
  //   double wei = res.first;
  //   vector<Lorentz5Momentum> p = res.second;
  //   PS += wei;
  //   vector<Energy>::const_iterator d = massvec.begin();
  //   if(1){
  //     cout << iev <<"   "<< wei << endl;
  //     for ( vector<Lorentz5Momentum>::iterator k = p.begin();
  //     k != p.end(); ++k, ++d ) {
  //       double px = k->x()/GeV; double py = k->y()/GeV; double pz = k->z()/GeV;
  //       double E = k->t()/GeV; double m = k->m()/GeV;
  //       cout << E <<"  "<< px <<"  "<< py <<"  "<< pz <<"  "<< m << endl;
  //     }
  //   }

  //   //cout << wei << endl;
  // }
  // cout << PS/nev << endl;

}




