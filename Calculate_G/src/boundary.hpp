#ifndef _BOUNDARY_HPP_
#define _BOUNDARY_HPP_

#include <dolfin.h>
using namespace dolfin;


class BoundaryX0 : public SubDomain {
    public:
    const double tol = 1e-14;
    double l,w,h;
    void setLwh( double lin, double win, double hin ) {
        l = lin; w = win; h = hin;
    }
    bool inside(const Array<double>& x, bool on_boundary) const override {
        return on_boundary && near(x[0],0,tol);
    }

};

class BoundaryX1 : public SubDomain {
    public:
    const double tol = 1e-14;
    double l,w,h;
    void setLwh( double lin, double win, double hin ) {
        l = lin; w = win; h = hin;
    }

    bool inside(const Array<double>& x, bool on_boundary) const override {
        return on_boundary && near(x[0],l,tol);
    }

};

class BoundaryY0 : public SubDomain {
    public:
    const double tol = 1e-14;
    double l,w,h;
    void setLwh( double lin, double win, double hin ) {
        l = lin; w = win; h = hin;
    }

    bool inside(const Array<double>& x, bool on_boundary) const override {
        return on_boundary && near(x[1],0,tol);
    }

};

class BoundaryY1 : public SubDomain {
    public:
    const double tol = 1e-14;
    double l,w,h;
    void setLwh( double lin, double win, double hin ) {
        l = lin; w = win; h = hin;
    }

    bool inside(const Array<double>& x, bool on_boundary) const override {
        return on_boundary && near(x[1],w,tol);
    }

};

class BoundaryZ0 : public SubDomain {
    public:
    const double tol = 1e-14;
    double l,w,h;
    void setLwh( double lin, double win, double hin ) {
        l = lin; w = win; h = hin;
    }

    bool inside(const Array<double>& x, bool on_boundary) const override {
        return on_boundary && near(x[2],h,tol);
    }

};

class BoundaryZ1 : public SubDomain {
    public:
    const double tol = 1e-14;
    double l,w,h;
    void setLwh( double lin, double win, double hin ) {
        l = lin; w = win; h = hin;
    }

    bool inside(const Array<double>& x, bool on_boundary) const override {
        return on_boundary && near(x[2],0,tol);
    }

};

#endif
