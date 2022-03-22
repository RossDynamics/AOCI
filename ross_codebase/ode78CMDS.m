function [tout, yout] = ode78CMDS(ode,tspan,y0,options,varargin)
%A wrapper for ode78 which makes it compatible with CMDS, tspan should be a
%two-element array in this case. If ode78 is in use, c.s.i.odeopts should
%set to a struct with fields called tol and trace.

[tout,yout] = ode78(ode, tspan(1), tspan(end), y0, tol, trace)