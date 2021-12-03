%   File: main.m
%   Copyright (c) <2020> <University of Mons>
%   COLORAMAP Group, Univ. of Mons, Belgium
%   https://sites.google.com/site/nicolasgillis/projects/overview
%
%   Permission is hereby granted, free of charge, to any person
%   obtaining a copy of this software and associated documentation
%   files (the "Software"), to deal in the Software without restriction,
%   including without limitation the rights to use, copy, modify and
%   merge the Software, subject to the following conditions:
%
%   1.) The Software is used for non-commercial research and
%       education purposes.
%
%   2.) The above copyright notice and this permission notice shall be
%       included in all copies or substantial portions of the Software.
%
%   3.) Publication, Distribution, Sublicensing, and/or Selling of
%       copies or parts of the Software requires special agreements
%       with the University of Mons and is in general not permitted.
%
%   4.) Modifications or contributions to the software must be
%       published under this license. The University of Mons
%       is granted the non-exclusive right to publish modifications
%       or contributions in future versions of the Software free of charge.
% 
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
%   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
%   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
%   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
%   OTHER DEALINGS IN THE SOFTWARE.
%
%   Persons using the Software are encouraged to notify the 
%   COLORAMAP Group, Univ. of Mons, Belgium about bugs. Please reference 
%   the Software in your publications if it was used for them.
%
% ------------------------------------------------------------------------
% SYNTAX:  
% This is a script not a function.  Parameters are set within. 
%
% OVERVIEW:
% This script runs all of the experiments and reproduces the plots from the
% associated paper.
%
% ------------------------------------------------------------------------
% REFERENCES:
% If this code is useful for you, please cite the paper:
% [1] 	Marrinan, Tim, P.-A. Absil, and Nicolas Gillis
%	"On a minimum enclosing ball of a collection of linear subspaces."
%     	arXiv preprint arXiv:2003.12455 (2020).
%
% ------------------------------------------------------------------------
% CREATED:      13/03/2020 by Tim Marrinan
%
% LAST EDITED:  13/03/2020 by Tim Marrinan
%
% NOTES: 
%
% 13/03/2020:
% The experiments in the paper are each run for 100 Monte Carlo trials, but
% to reproduce this will take hours, like maybe 20(?) hours, for all plots.
% Just be warned.  You can see that the behavior is consistent from a much
% smaller number of trials.
% ------------------------------------------------------------------------

scen{1}     = 'Fig. 3a';
scen{2}     = 'Fig. 3b';
scen{3}     = 'Fig. 4a';
scen{4}     = 'Fig. 4b';
scen{5}     = 'Fig. 5';
scen{6}     = 'Fig. 6';
scen{7}     = 'Fig. 7a';
scen{8}     = 'Fig. 7b';
scen{9}     = 'Fig. 8';
scen{10}    = 'Fig. 9';
scen{11}    = 'Fig. 10';
nRuns       = 5;
queued      = 1:11;

for i = queued
    scenario = scen{i};
    if i <= 4
        GMEB_IllustrativeExample(scenario);
    elseif i > 4 && i <= 6
        GMEB_Exp001(scenario,nRuns);
    elseif i > 6 && i <=  8
        GMEB_Exp002(scenario,nRuns);
    else
        GMEB_Exp003(scenario,nRuns);
    end
end
fprintf('\tFinished.\n');






