#ifndef _GKRONROD_HPP_
#define _GKRONROD_HPP_
#include <valarray>
#include <array>
#include <set>
#include <algorithm>

using namespace std;

template <typename val_type,typename fn_type>
valarray<val_type> myapply(valarray<val_type> xs,fn_type f);




template <typename real_type>
class GKData{
  private:
    valarray<real_type> _gknodes = {-0.991455371120813,-0.949107912342759,-0.864864423359769,-0.741531185599394,-0.586087235467691,-0.405845151377397,-0.207784955007898,0.0,0.207784955007898,0.405845151377397,0.586087235467691,0.741531185599394,0.864864423359769,0.949107912342759,0.991455371120813};
    valarray<real_type> _kweights = {0.022935322010529,0.063092092629979,0.104790010322250,0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728,0.204432940075298,0.190350578064785,0.169004726639267,0.140653259715525,0.104790010322250,0.063092092629979,0.022935322010529};
    valarray<real_type> _gweights = {0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469};

  public:
    valarray<real_type> gknodes() {return  this->_gknodes ;};
    valarray<real_type> kweights(){return  this->_kweights;};
    valarray<real_type> gweights(){return  this->_gweights;};

};

template <typename real_type, typename fn_type>
class interval{
  private:
    real_type quadval;
    real_type errval;
  public:
    array<real_type,2> endpoints;
    interval(real_type left, real_type right){
      this->quadval=0.0;
      this->errval =1.0;
      this->endpoints[0]=left;
      this->endpoints[1]=right;
    }
    void gkquad(fn_type f,valarray<real_type> nodes, valarray<real_type> kweights,valarray<real_type> gweights){
      /*Map nodes from [-1,1] to [endpoints[0],endpoints[1]].*/
      auto a        = endpoints[0];
      auto b        = endpoints[1];
      auto nds      = myapply(nodes,[a,b](real_type x)->{return (a-x)*0.5+(b+x)*0.5;});
      auto h        = b-a;
      /*Scale weights.*/
      auto wts      = kweights*h/2.0;
      auto gwts     = gweights*h/2.0;
      /*Evaluate function at nodes.*/
      auto fx       = myapply<real_type,fn_type>(nds,f);
      /*Compute Kronrod quadrature.*/
      this->quadval = (fx*wts).sum();

      /*Calculate Gauss quadrature.*/
      slice slc(1,15,2);
      valarray<real_type> slcfx = fx[slc];
      int n=nds.size();
      for(int i=0;i<n;i++)
        cout <<endl<<nds[i]<<endl;
      auto gquadval = (slcfx*gwts).sum();
      /*This error estimate is a little suspicious..*/
      this->errval = 200.0*pow(abs(quadval-gquadval),1.5);
    }
    real_type get_error(){
      return this->errval;
    }
    real_type get_quad(){
      return this->quadval;
    }
    real_type left(){return this->endpoints[0];}
    real_type right(){return this->endpoints[1];}



};


template <typename real_type,typename val_type>
val_type gauss_kronrod(val_type (*f)(real_type),real_type reltol=1e-3);


template <typename real_type, typename fn_type>
bool interval_lessthan(const interval<real_type,fn_type>& int1, const interval<real_type,fn_type>& int2);



template <typename real_type,typename fn_type>
real_type gauss_kronrod(fn_type f,real_type reltol){
  /*Get gauss-kronrod data.*/
  GKData<real_type> dat;
  auto gknodes  = dat.gknodes ();
  auto kweights = dat.kweights();
  auto gweights = dat.gweights();

  /*Start algorithm at reference interval [-1,1].*/
  set<interval<real_type, fn_type>, bool(*)(const interval<real_type,fn_type>&,const interval<real_type,fn_type>&) > intervals(&interval_lessthan<real_type,fn_type>);
  intervals.insert(interval<real_type,fn_type>(-1.0,1.0));
  bool has_split=true;
  while (has_split){
    has_split=false;
    auto tmp_intervals = intervals;
    for(auto ab : tmp_intervals){
      ab.gkquad(f,gknodes,kweights,gweights);
      if(ab.get_error()>reltol){
        has_split=true;
        auto lp = ab.left(); 
        auto rp = ab.right();
        auto mp = (lp+rp)/2.0;
        intervals.insert(interval<real_type,fn_type>(lp,mp));
        intervals.insert(interval<real_type,fn_type>(mp,rp));
        intervals.erase(ab);
      }
    }
  }


  size_t nints = intervals.size();
  size_t intid = 0;
  valarray<real_type> quadvals(nints);
  for(auto ab : intervals){
    quadvals[intid]=ab.get_quad();
    intid++;
  }

  return quadvals.sum();
}

template <typename val_type,typename fn_type>
valarray<val_type> myapply(valarray<val_type> xs,fn_type f){
  size_t len = xs.size();
  valarray<val_type> fxs(len);
  for(size_t it=0;it<len;it++)
    fxs[it]=f(xs[it]);
  return fxs;
}

template <typename real_type,typename fn_type>
bool interval_lessthan(const interval<real_type,fn_type>& int1, const interval<real_type,fn_type>& int2){
  return lexicographical_compare(int1.endpoints.begin(),
      int1.endpoints.end(),
      int2.endpoints.begin(),
      int2.endpoints.end());
}



#endif
