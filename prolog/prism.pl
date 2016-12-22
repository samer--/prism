:- module(prism, [
  % session management
		prism_start/0,	prism_start/1,	prism_start/2
	,	prism_restart/0
	,	prism_recover/0
	,	prism_close/0
   ,  prism_is_running/0

	% state
	,	prism_state_get/1
	,	prism_state_set/1

	% loading source files
	,	load_bprolog/1
	,	load_prism/1,	load_prism/2

	% sampling
	,	prism_sample/1
	,	prism_sample_n/3

	% inference
	,	prob/2
	,	explain/2
   ,  explain_probs/3
	,	hindsight/3
	,	hindsight_agg/3
	,	chindsight/3
	,	chindsight_agg/3
	,	viterbi/4, viterbi_expl/5, viterbi_tree/5
	,	viterbi_graph_to_tree/2

	% learning
	,	prism_learn/4

	% pretty printing
	,	print_graph/1
	,	print_tree/1

	% switch management
	,	sw_set/3
	,	sw_get/3 
	,	sw_fix/2
	,	sw_unfix/2
	,	sw_values/2
	,	sw_set_to_sample/1, sw_set_to_sample/2
	,	sw_set_to_mean/1, sw_set_to_mean/2

	% flags
	,	prism_flag/3
	,	prism_flag_affects/2
	,	prism_flag_set/2
	,	prism_flag_get/2

	% information
	,	prism_statistics/3,	prism_statistics/2
	,	prism_show/1

	% low level interface
	,	prism/1
	,	prism_nd/1
	,	(prism_dynamic)/1

	,	(#)/1
	,	(##)/1

	% operators
	,	op(1150,fx,prism_dynamic)
	,	op(950,fx,#)
	,	op(950,fx,##)
	]).

/** <module> Using PRISM as a child process

This module provides access PRISM, a Prolog-based probabilistic programming
system. PRISM runs under B-Prolog, so this module provides ways to manage
and communicate with a PRISM/B-Prolog child process.

Much of this can be seen as an effort to make the interface less stateful.
Currently, the state of the probabilistic switches (distributions and pseudocounts) 
is the main piece of state. Most of the flags that affect learning and inference 
are managed explicitly and statelessly by the relevant procedures.

Flags which are still effective statefully are:
	- flags to do with information and progress display
	- default parameters for switches
	- explanation graph housekeeping
	- sort_hindsight
	- log_scale affects explanation, viterbi and learning

Types used in this module:
==
prism_goal   % a callable term
prism_state  % term describing state of PRISM 
filename     % a literal absolute or relative filename 
filepath     % path to file using search path mechanism.
==

*/
:- use_module(library(plrand), [with_rnd_state/1]).

:- dynamic current_prism/3.
:- dynamic current_prism_file/1.
:- dynamic (prism_dynamic)/2.
:- dynamic saved_state/2.

:- prolog_load_context(directory, Dir),
   directory_file_path(Dir,'psm',PSM),
   assert(user:file_search_path(prism,PSM)).



% ------------------------------------------------------------------------
%                             MAIN INTERFACE 
% ------------------------------------------------------------------------

% ---------------------- SESSION CONTROL----------------------------------

%% prism_start(+Exec:atom, +LogFile:filename) is det.
%% prism_start(+LogFile:filename) is det.
%% prism_start is det.
%
%  Start PRISM as a child process. Exec should be a standard filespec pointing to the 
%  PRISM executable.  %  Output from PRISM is recorded to LogFile. PRISM immediately 
%  runs bprepl.pl, which receives queries on stdin and replies on stdout. Any currently 
%  running PRISM %  is closed first. If not supplied, executable is set to =|path(prism)|=.
%  If not supplied, log file is 'prism.log'. 
prism_start :- prism_start('prism.log').
prism_start(Log) :- 
   prism_start(path(prism), Log).
prism_start(ExecSpec,Log) :-
	absolute_file_name(ExecSpec, [access(execute)],Exec),
	absolute_file_name(prism('bprepl.pl'),[access(read)],S),
	format(atom(Args),"cl('~w'),main(file('~w'))",[S,Log]),
   start_prism_with_args(Exec, ["-g", Args]).

start_prism_with_args(Exec,Args) :-
	process_create(Exec,Args,[stdin(pipe(In)),stdout(pipe(Out)),process(PID)]),
	(current_prism(_,_,_) -> prism_close; true),
	set_streams(In,Out,child(Exec,Args,PID)),
	nl, once((repeat, wait(infinite), recv(term(repl(ready),_)))),
	foreach(prism_dynamic(FF,AA),prism(dynamic(FF/AA))).

set_streams(In,Out,Descriptor) :-
	%	set_stream(Out,timeout(0)),
	set_stream(Out,close_on_abort(false)),
	set_stream(In,close_on_abort(false)),
	assert(current_prism(Descriptor,In,Out)),
	set_prolog_flag(float_format,'%.15g'). % for accurate exchange of floats


%% prism_is_running is semidet.
%
%  Succeeds if an instance of prism is running. 
prism_is_running :- current_prism(_,_,_).

%% prism_restart is det.
%
%  Restart the current PRISM process, which must be a child as created
%  by prism_start/0 or prism_start/2. The current state is saved and
%  restored to the new PRISM process.
prism_restart :-
	prism_full_state(S),
	prism_start_with(S).

%% prism_restore_state(+ID) is det.
%  Retrieves PRISM state previously stored under ID and
%  installs it into PRISM. See prism_save_state/1.
prism_restore_state(ID) :-
	saved_state(ID,S),
	prism_start_with(S).


%% prism_save_state(+ID) is det.
%  Get current PRISM state and save it under the name ID.
%  Can later be restored used prism_restore_state(ID).
prism_save_state(ID) :-
	prism_full_state(S),
	retractall(saved_state(ID,_)),
	assert(saved_state(ID,S)).

prism_full_state(state(Exec,Args,Files,State)) :-
	prism_state_get(State),
	findall(F,current_prism_file(F),Files),
	current_prism(child(Exec,Args,_),_,_).


%% prism_recover is det.
%
%  Try to restart from saved state if there is one.
%  Recovery state can be saved from prism_restart/0 and
%  prism_restore_state/1. 

prism_recover :- 
	recovery_state(nothing,nothing), !,
	format('Nothing to recover from\n').

prism_recover :-
	recovery_state(just(State),nothing),
	prism_start_with(State).

recovery_state(S1,S2) :-
	( saved_state(recovery,S) -> S1=just(S); S1=nothing),
	retractall(saved_state(recovery,_)),
	(	S2=just(SS2) 
	->	assert(saved_state(recovery,SS2))
	;	true).


%% prism_start_with(+S:prism_state) is det.
%
%  Tries to start PRISM with the given state, saving
%  the state as recovery_state if it fails.
prism_start_with(S) :-
	S=state(Exec,Args,Files,State),
	catch((
			start_prism_with_args(Exec,Args),
			maplist(prism,Files),
			prism_state_set(State)
		), Ex, (
			recovery_state(_,just(S)),
			format('PRISM exception: ~w.\n',[Ex]),
			format('State was saved. Use prism_recover when the problem has been fixed.\n'),
			throw(Ex)
		)
	).


%% prism_close is det.
%
%  Close the current PRISM process.
prism_close :-
	(	current_prism(Desc,In,Out) 
	-> writeln('Closing PRISM'),
		catch((close(In), close(Out)), Ex, print_message(informational,Ex)),
      (Desc=child(_,_,PID) -> process_wait(PID,_); true),
		retractall(current_prism(_,_,_))
	;	writeln('No active PRISM.')
	).

% ---------------------- Low level interface ----------------------------
 
%% prism( Goal:prism_goal) is semidet.
%
%  Submit query to PRISM process and get reply. Goal is a normal Prolog goal which
%  submitted to PRISM. It is run as by once/1 and any bindings returned.
prism(Goal) :- 
	gensym(q,Id),
	term_variables(Goal,Vars),
	send(query(Id,Goal,Vars)), wait(infinite), 
	once((repeat,recv(term(query(Id,Vars,Reply),_)))),
	interp(Reply).

%% prism_nd( Goal:prism_goal) is nondet.
%
%  Submit nondeterministic query to PRISM. The goal in run in PRISM via
%  findall/3. All the results are returned and retrieved on backtracking
%  in this predicate.
prism_nd(Q) :- term_variables(Q,V), prism(findall(V,Q,VX)), member(V,VX).


%% #(Goal:prism_goal) is semidet.
%  Unary prefix operator, equivalent to prism/1.
#(G) :- prism(G).

%% ##(Goal:prism_goal) is semidet.
%  Unary prefix operator, equivalent to prism_nd/1.
##(G) :- prism_nd(G).



% ------------------------------------------------------------------------
%                   Higher level prism operations 
% ------------------------------------------------------------------------

% ---------------------- PRISM state -------------------------------------

%% prism_state( -S1:prism_state, +S2:prism_state) is det.
%
%  Perform a transition of the global PRISM state, retreiving the current
%  state into S1 and then setting it to S2.
prism_state(S1,S2) :- prism_state_get(S1), prism_state_set(S2).

%% prism_state_get( -S:prism_state) is det.
%  Get the current state of the PRISM system into a Prolog term S. Can
%  be restored using prism_state_set/1.
prism_state_get(ps(SX,PX,CX,FX,VX,DX)) :-
	prism(findall(sw(I,V), (get_reg_sw(I), get_values1(I,V)),SX)),           % get all switch values
	prism(findall(sw(I,V), dynamic_values(I,V),VX)),   % just the dynamic ones
	prism(findall(sw(I,p,set(F,P)),(get_reg_sw(I),get_sw(I,[F,_,P])),PX)),   % all switch probabilities
   % !! NB I changed this from h to a as h is a synonym for d, not a
	prism(findall(sw(I,a,set(F,C)),(get_reg_sw(I),get_sw_a(I,[F,_,C])),CX)), % all switch pseudocounts
	prism(findall(flag(F,V), get_prism_flag(F,V),FX)), % all prism flags
	findall( dyn(FF,AA,CCX), (
			prism_dynamic(FF,AA), functor(PP,FF,AA),
			prism(findall(PP,PP,CCX))
		), DX).
		

%% prism_state_set( +S:prism_state) is det.
%  Set the current state of the PRISM system to that described in S. State
%  can be obtained using prism_state_get/1 or prism_state/2.
prism_state_set(ps(SX,PX,CX,FX,VX,DX)) :-
	foreach((member(flag(F,V),FX), V\='$disabled'), prism_flag_set(F,V)),       % restore flags
	foreach(member(sw(I,V),VX), prism(set_dynamic_values(I,V))), % assert dynamic values if necessary
	foreach(member(sw(I,V),SX), check(prism(get_values(I,V)))),      % check all switch values
	forall( member(dyn(_,_,CCX),DX), 
		prism(forall(member(CC,CCX),assert(CC)))),
	maplist(set_switch,PX), % set switch probabilities
	maplist(set_switch,CX). % set switch pseudocounts

sw_info(probs,sw(I,_),sw(I,p,P)) :- prism_nd(get_sw(I,[F,_,PX])) *-> P=set(F,PX); P=unset.
sw_info(counts,sw(I,_),sw(I,a,C)) :- prism_nd(get_sw_a(I,[F,_,CX])) *-> C=set(F,CX); C=unset.
sw_info(exponents,sw(I,_),sw(I,h,C)) :- prism_nd(get_sw_d(I,[F,_,CX])) *-> C=set(F,CX); C=unset.

set_switch(sw(_,_,unset)) :- !. % no way to explicity un-set parameters
set_switch(sw(I,p,set(F,P))) :- !, prism(get_sw(I,P)), (F=fixed->prism(fix_sw(I));prism(unfix_sw(I))).
set_switch(sw(I,a,set(F,C))) :- !, 
   prism(set_sw_a(I,C)), (F=fixed_a->prism(fix_sw_a(I));prism(unfix_sw_a(I))).
set_switch(sw(I,h,set(F,C))) :- !, 
   maplist(add(1),C,A), prism(set_sw_a(I,A)), % use set_sw_a incase pseudocounts are negative
   (F=fixed_a->prism(fix_sw_a(I));prism(unfix_sw_a(I))).
set_switch(I) :- format('*** unrecognised switch setting: ~q.\n',[I]).

check(G) :- call(G) -> true; format('*** check failed: ~w.\n',[G]), fail.

%% prism_dynamic(+Predicates) is det.
%
%  Registers predicates as dynamic predicates in PRISM. These
%  are then declared as dynamic automatically whenever PRISM is started.
prism_dynamic(F/A) :- prism_dynamic(F,A), !.
prism_dynamic(F/A) :- !, assert(prism_dynamic(F,A)).
prism_dynamic((P,PX)) :- prism_dynamic(P), prism_dynamic(PX).

stoch(H,P,N) :- sumlist(H,N), N>0, maplist(divby(N),H,P).
divby(N,X,Y) :- Y is X/N.


% ---------------------- PROGRAM LOADING ----------------------------

%% load_bprolog( +Name:filename) is det. 
%  Compile and link a B-Prolog source file.
load_bprolog(Name)  :- prism(cl(Name)), assert(current_prism_file(cl(Name))).

%% load_prism( +Name:filepath, +Opts:list(oneof([neg]))) is det.
%% load_prism( +Name:filepath) is det.
%
%  Compile and link a PRISM source file. If file extension is 
%  ommitted, then '.psm' is assumed. There is only one option currently:
%  if neg is included in the option list, the program is assumed to
%  include negation a failure predicate, and is loaded via prismn/1.
load_prism(Name) :- load_prism(Name, []). 
load_prism(Name,Opts) :- 
	absolute_file_name(Name,[access(read),extensions([psm])],Filename),
	file_directory_name(Filename,Directory),
	prism(cd(Directory)),
	(  member(neg,Opts) 
	-> prism(prismn(Filename))
	;  prism(prism(Filename))
	),
	retractall(current_prism_file(prism(_))), % only one prism file active at a time
	assert(current_prism_file(prism(Name))).

% ---------------------- MANAGING SWITCHES ----------------------------

%% sw_values(+S:switch,-V:list(A)) is semidet.
%% sw_values(-S:switch,-V:list(A)) is nondet.
%  True when S is a registered switch and V is its list of values.
sw_values(S,V)    :- prism_nd((get_reg_sw(S),get_values1(S,V))).

%% sw_get(+Spec:oneof([probs,counts]), +S:switch, -I:switch_info) is det.
%% sw_get(-Spec:oneof([probs,counts]), -S:switch, -I:switch_info) is nondet.
%
%  Gets information on probability distribution or pseudocounts for
%  a switch. The returned switch_info type is
%  ==
%  switch_info ---> unset ; set(oneof([fixed,unfixed]), list(float)).
%  ==
sw_get(T,S,I)  :- sw_values(S,V), sw_info(T,sw(S,V),sw(_,_,I)).

%% sw_set(+Spec, +S:switch, +H:list(float)) is det.
%
%  Set a switch's parameters depending on Spec, which can be:
%  * dist 
%    H is an unnormalised distribution. It will be normalised.
%  * probs 
%    H is a normalised distribution, and will be installed.
%  * counts 
%    H will be set as the switch's pseudocounts.
%  * all_probs
%    Sets probabilities for all switchs that unify with S
%  * all_counts
%    Sets counts for all switchs that unify with S
sw_set(dist,S,H)  :- stoch(H,V,_), prism(get_sw(S,V)).
sw_set(probs,S,V)  :- prism(get_sw(S,V)).
sw_set(counts,S,V) :- prism(set_sw_a(S,V)).
sw_set(all_probs,S,V) :- prism(set_sw_all(S,V)).
sw_set(all_counts,S,V) :- prism(set_sw_all_a(S,V)).

%% sw_fix(+Spec:oneof([prob,counts]), S:switch) is det.
%  Fixes the distribution or pseudocounts for the named switch.
%  This prevents the setting from changing during learning.
sw_fix(probs,S)    :- prism(fix_sw(S)).
sw_fix(counts,S)    :- prism(fix_sw_a(S)).

%% sw_unfix(+Spec:oneof([prob,counts]), S:switch) is det.
%  Unfixes the distribution or pseudocounts for the named switch.
%  This allows the setting to change during learning.
sw_unfix(probs,S)  :- prism(unfix_sw(S)).
sw_unfix(counts,S)  :- prism(unfix_sw_a(S)).


%% switch_distribution(-Type) is nondet.
%  Table of distribution types for switch setting. 
switch_distribution(default).
switch_distribution(uniform).
switch_distribution(f_geometric).
switch_distribution(f_geometric(natural)).
switch_distribution(f_geometric(natural,oneof([asc,desc]))).

%% sw_set_to_mean(+S:switch) is det.
%% sw_set_to_mean(+S:switch,+Offset:float) is det.
%
%  Sets a switch's distribution to the mean of the
%  Dirichlet distribution whose parameters come from the
%  current 'counts' setting of the switch. Off, if present, is 
%  added to all the pseudocounts.
sw_set_to_mean(S) :- sw_set_to_mean(S,0).
sw_set_to_mean(S,Off) :-
	sw_get(counts,S,set(_,D)),
	maplist(add(Off+1),D,A),
	sw_set(dist,S,A).

%% sw_set_to_sample(+S:switch) is det.
%% sw_set_to_sample(+S:switch, Offset:float) is det.
%
%  Sample a new discrete distribution for a switch's parameters
%  from a Dirichlet distribution whose parameters come from the
%  current 'counts' setting of the switch. Off, if present, is 
%  added to all the pseudocounts. This uses the undocumented 
%  sample_Dirichlet/5 predicate from library(plrand).
sw_set_to_sample(S) :- sw_set_to_sample(S,0).
sw_set_to_sample(S,Off) :-
	sw_get(counts,S,set(_,D)),
	maplist(add(Off+1),D,A),
   length(A,N),
	with_rnd_state(plrand:sample_Dirichlet(N,A,P)),
	sw_set(probs,S,P).

add(Z,X,Y) :- Y is X+Z. 


% ---------------------- SAMPLING EXECUTION ----------------------------

%% prism_sample(Goal:prism_goal) is semidet.
%  Sample from distribution specified by PRISM goal.
prism_sample(Q) :- prism(sample(Q)).

%% prism_sample_n(+N:nat, Goal:cond_goal, Results:list(prism_goal)) is det.
%
%  Multiple samples from distribution specified by PRISM goal.
%  Goal can be an ordinary PRISM goal or a conditional goal of the
%  form (Goal|Cond), where Cond is an ordinary Prolog goal.
%  Returns a list of copies of the called goal with return values
%  filled in.
prism_sample_n(N,Q,X) :- Q\=(_|_), !, prism(get_samples(N,Q,X)).
prism_sample_n(N,(Q|C),X) :- prism(get_samples_c(N,Q,C,X)).

% ---------------------- INFERENCE ----------------------------

%% prob(Goal:prism_goal, -P:float) is det.
% Compute probability of given goal.
prob(Q,P) :- prism(prob(Q,P)).

%% explain(Goal:prism_goal,-Graph) is det.
%  Gets explanation graph for given goal.
explain(Q,G)          :- prism(probf(Q,G)).

%% explain_probs(+Type:oneof([inside,outside,viterbi,in_out]), Goal:prism_goal,-Graph) is det.
%  Gets explanation graph annoted with probabilities for given goal. The explanation
%  types map to a PRISM commands as follows:
%     * inside -> probfi/2
%     * outside -> probfo/2
%     * in_out -> probfio/2
%     * viterbi -> probvf/2
explain_probs(inside,Q,G)   :- prism(probfi(Q,G)).
explain_probs(outside,Q,G)  :- prism(probfo(Q,G)).
explain_probs(viterbi,Q,G)  :- prism(probfv(Q,G)).
explain_probs(in_out,Q,G)   :- prism(probfio(Q,G)).

%% hindsight(TopGoal:prism_goal, SubGoal:prism_goal,-Probs:list(pair(prism_goal,nonneg))) is nondet.
%
% Computes probabilities of subgoals unifying with SubGoal for a given top goal. Matches
% different subgoals on backtracking.
% hindsight(Q,R,P)      :- prism(hindsight(Q,R,P1)), maplist(list_pair,P1,P).
hindsight(Q,R,P) :-      prism(hindsight(Q,R,P1)), member([R,P],P1).
chindsight(Q,R,P) :-     prism(chindsight(Q,R,P1)), member([R,P],P1).

%% hindsight_agg(TopGoal:prism_goal, SubGoal:prism_goal,-Probs:list(list(pair(prism_goal,nonneg)))) is det.
%
% Computes total probability of all subgoals unifying with SubGoal for a given top goal.
% NB. this functionality can be replicated in a more idiomatic way using hindsight/3 and
% SWI Prolog's library(aggregate).
hindsight_agg(Q,R,P) :- prism(hindsight_agg(Q,R,P1)), maplist((maplist(list_pair)),P1,P).

%% chindsight(TopGoal:prism_goal, SubGoal:prism_goal,-Probs:list(pair(prism_goal,nonneg))) is nondet.
%
% Computes conditional probabilities of matching subgoals given that the top goal is true.
% Same as hindsight/3 but with probabilities divided by the probabality of the top goal.
% chindsight(Q,R,P)     :- prism(chindsight(Q,R,P1)), maplist(list_pair,P1,P).

%% chindsight_agg(TopGoal:prism_goal, SubGoal:prism_goal,-Probs:list(list(pair(prism_goal,nonneg)))) is det.
%
% Computes aggregate conditional probability of matching subgoals given that the top goal is true.
% Same as chindsight_agg/3 but with probabilities divided by the probabality of the top goal.
% NB. this functionality can be replicated in a more idiomatic way using hindsight/3 and
% SWI Prolog's library(aggregate).
chindsight_agg(Q,R,P) :- prism(chindsight_agg(Q,R,P1)), maplist((maplist(list_pair)),P1,P).

goal_probs(P1,P2) :- maplist(list_pair,P1,P2).
list_pair([X,Y],X-Y).

% expectation
expect(X,G,EX) :-
	chindsight(G,G,PX),
	accum(X,G,PX,EX).

accum(_,_,[], 0).
accum(X,G,[G1-P1 | PX], Tot) :-
	accum(X,G,PX,TotX),
	copy_term(X/G,X1/G1),
	Tot is P1*X1 + TotX.


%% viterbi(+N:nat, Goal:prism_goal, -P:float, +Opts:list(viterbi_opt)) is nondet.
%  Compute probabilities of N most probable explanation for Goal.
%  Options can be as follows:
%  ==
%  viterbi_opt ---> ground(bool)        % [false] whether or not to ground variable in Goal
%                 ; mode(oneof([ml,vb]) % [ml]  use probs (ml) or counts (vb) 
%                 ; rerank(nat)         % [10]  number of explanations for reranking
%                 .
%  ==

viterbi(N,Q,P,Opts) :- 
	process_viterbi_opts(Opts,Ground),
	viterbi_prob(Ground,N,Q,P).

viterbi_prob(false,1,Q,P) :- !, prism(viterbi(Q,P)).
viterbi_prob(false,N,Q,P) :- prism(n_viterbi(N,Q,PX)), member(P,PX).
viterbi_prob(true,1,Q,P)  :- !, prism(viterbig(Q,P)).
viterbi_prob(true,N,Q,P)  :- prism(n_viterbig(N,Q,P)).


%% viterbi_expl(+N:nat, Goal:prism_goal, -P:float, +Opts:list(viterbi_opt), -G:graph) is nondet.
%  Compute probabilities of N most probable explanation for Goal and
%  also produce graphs for each. See viterbi/4 for options.
viterbi_expl(N,Q,P,Opts,G) :- 
	process_viterbi_opts(Opts,Ground),
	viterbi_expl1(Ground,N,Q,P,G).

viterbi_expl1(false,1,Q,P,G) :- !, prism(viterbif(Q,P,G)).
viterbi_expl1(false,N,Q,P,G) :- prism(n_viterbif(N,Q,VX)), member(v_expl(_,P,G),VX).
viterbi_expl1(true,1,Q,P,G)  :- !, prism(viterbig(Q,P,G)).
viterbi_expl1(true,N,Q,P,G)  :- prism(n_viterbig(N,Q,P,G)).

%% viterbi_tree(+N:nat, Goal:prism_goal, -P:float, +Opts:list(viterbi_opt), -T:tree) is nondet.
%  Compute probabilities of N most probable explanation for Goal, as viterbi_expl/5,
%  but produce a tree instead of a graph by using viterbi_graph_to_tree/2.
%  also produce tree for each. See viterbi/4 for options.
viterbi_tree(N,Q,P,Opts,T) :- 
   viterbi_expl(N,Q,P,Opts,G),
   viterbi_graph_to_tree(G,T).

process_viterbi_opts(Opts,Ground) :-
	select_option(ground(Ground),Opts,O1,false),
	select_option(mode(Mode),O1,O2,ml),
	select_option(rerank(Rerank),O2,O3,10),
	(	O3=[] -> true; throw(error(unrecognised_options(O2)))),
	prism_flag_set(viterbi_mode,Mode),
	prism_flag_set(rerank,Rerank).


%% viterbi_graph_to_tree( +G:graph, -T:tree) is det.
%
%  Computes Viterbi tree for a given Viterbi explanation graph.
%  Graph comes from one of the viterbi or n_viterbi predicates.
viterbi_graph_to_tree(G,T) :- prism(viterbi_tree(G,T)).

%% print_graph(+G:graph) is det.
%  Prints an explanation graph on PRISM's user output stream.
print_graph(G) :- prism(print_graph(user_output,G,[])).

%% print_tree(+G:tree) is det.
%  Prints a Viterbi tree PRISM's user output stream.
print_tree(T)  :- prism(print_tree(user_output,T,[])).



% ---------------------- LEARNING ------------------------

%% prism_learn( +Method:learn_method, +GG:list(goal), +Opts:list(learn_opt), -Scores:list(learn_score)) is det.
%
% Learn model parameters from a list of goals.
%
% ==
%  learn_method ---> map(Init:oneof([none,random,noisy(Std)]))
%                  ; vb(Init:oneof([none,perturb(Std),reset,noisy(Std),Pt:oneof([on,off]))
%                  ; vb_pm(Init:oneof([none,perturb(Std),reset,noisy(Std))
%                  .
% ==
% For MAP learning, the Init option determines how switch probabilities are
% initialised; the values have the following meanings:
%   * none
%     Switch probabilities begin with their current values.
%   * random
%     Probabilities are set to random values
%   * noisy(Std:nonneg)
%     Probabilities are drawn from rectified Gaussian with given standard deviation.
%
% For VB learning, initialisation methods are as follows
%
%   * none
%     Learning continues from current pseudocounts (ie evidence accumulates).
%   * perturb(Std:nonneg)
%     Current values have rectified Gaussian noise added
%   * reset
%     Hyperparameters are set to value of default_sw_a flag.
%   * noisy(Std:nonneg)
%     Hyperparameters are set to value of default_sw_a plus some noise.
%
% Method vb_pm is like vb except that switch parameters are set to their posterior
% means after learning.
%
% Valid options are (with defaults in comments): 
% ==
%  learn_opt ---> daem(oneof([off,sched(Init,Rate)])) % on
%               ; epsilon(nonneg)                     % 0.0001
%               ; max_iterate(nonneg)                 % default
%               ; restart(natural)                    % 1
%               .
%
%  learn_score ---> log_prior(float)    % only for map learning
%                 ; log_post(float)     % only for map learning
%                 ; log_lik(float)
%                 ; free_energy(float)  % only for VB learning
%                 .
% ==

prism_learn(map(Init),GG,Opts,[log_prior(LPr),log_lik(LL),log_post(LPo),bic(BIC)]) :- 
	set_prism_option(learn_mode,ml),
	set_prism_option(init,Init),
	process_learn_options(Opts),
	prism_flush(learn(GG)),
   prism_statistics(log_likelihood,LL),
   (  maplist(prism_statistics,[log_prior,log_post,bic],[LPr,LPo,BIC]) -> true
   ;  prism_statistics(num_parameters,K),
      length(GG,N),
      LPr=0, LPo=LL, 
      BIC is LL - (K/2)*log(N)
   ).

prism_learn(vb(Init),GG,Opts,[free_energy(FE)]) :- 
	set_prism_option(vb_init,Init),
	set_prism_option(learn_mode,vb),
	process_learn_options(Opts),
	prism_flush(learn(GG)),
	prism_statistics(free_energy,FE).


prism_learn(vb_pm(Init),GG,Opts,[free_energy(FE)]) :-
	set_prism_option(vb_init,Init),
	set_prism_option(learn_mode,both),
	process_learn_options(Opts),
	prism_flush(learn(GG)),
	prism_statistics(free_energy,FE).


process_learn_options(Options) :-
	foldl(process_option, 
		[	fix_init_order	/on
		,	epsilon  		    /0.0001
		,	max_iterate 	  /default
		, daem            /off
		,	restart 			  /1
		],Options,O2),
		(	O2=[] -> true; throw(error(unrecognised_options(O2)))).

process_option(Name/Default,O1,O2) :- 
	functor(Opt,Name,1), arg(1,Opt,Val),
	select_option(Opt,O1,O2,Default),
	set_prism_option(Name,Val).


set_prism_option(vb_init,Init) :- !,
	member(Init/(Reset,Std),[none/(off,0.000000001), perturb(U)/(off,U), 
									 reset/(on,0.000000001), noisy(U)/(on,U)]),
	prism_flag_set(reset_hparams,Reset),
	prism_flag_set(std_ratio,Std).

set_prism_option(init,noisy(Std)) :- !,
	prism_flag_set(init,noisy_u), 
	prism_flag_set(std_ratio,Std).

set_prism_option(daem,off) :- !, prism_flag_set(daem,off).
set_prism_option(daem,sched(Bstart,Brate)) :- !, 
	prism_flag_set(daem,on),
	prism_flag_set(itemp_init,Bstart),
	prism_flag_set(itemp_rate,Brate).
	
set_prism_option(Name,Val) :- prism_flag_set(Name,Val).

% --------------------- PRISM FLAGS ---------------------------


%% prism_flag_set( +Name, +Value) is det.
%  Set value of the named PRISM flag.
prism_flag_set(F,V) :- prism(set_prism_flag(F,V)).

%% prism_flag_get( +Name, -Value) is det.
%% prism_flag_get( -Name, -Value) is nondet.
%
%  Get value of the named PRISM flag or of all flags in turn.
prism_flag_get(F,V) :- prism_nd(get_prism_flag(F,V)).

%% prism_flag( +Name, -Type, -Description) is det.
%% prism_flag( ?Name, ?Type, ?Description) is nondet.
%
%  Contains information about all the PRISM flags.

% cosmetic - affect display only
prism_flag( warn,       	oneof([on,off]), 'controls display of warning messages').
prism_flag( verb,      		   oneof([none,graph,em,full]), 'control verbosity in learning').
prism_flag( show_itemp,       oneof([on,off]), 'print inverse temperature rate in DAEM algorithm').
prism_flag( search_progress, 	natural,         'frequency of message printing in explanation finding').
prism_flag( em_progress,      natural,         'frequence of message printing in EM algorirhm').
prism_flag( mcmc_progress, 	natural,         'frequency of message printing in MCMC learning').
prism_flag( learn_message, 	oneof([search,em,stats,misc,none,all]),      'messages during learning (can combine with +)').
prism_flag( mcmc_message,   	oneof([search,em,mcmc,stats,misc,none,all]), 'messages during MCMC learning (can combine with +)').
prism_flag( write_call_events,  oneof([none, off, all, _]), 'no idea').

% model initialisation
prism_flag( default_sw_a,  oneof([none,uniform,uniform(nonneg)]), 'default values for switch pseudo-counts').
prism_flag( default_sw_d,  oneof([none,uniform,uniform(nonneg)]), 'default values for switch pseudo-counts').
prism_flag( default_sw, 	
	[none, uniform, f_geometric, f_geometric(natural), f_geometric(natural,oneof([asc,desc]))], 
	'default distribution for new switches').

% affect inference, explanation generation
prism_flag( clean_table, 	oneof([on,off]), 'dispose of explanation table after use').
prism_flag( force_gc,   	oneof([on,off]), 'trigger garbage collection after inference operations').
prism_flag( reduce_copy, 	oneof([on,off]), 'affects copying of terms in explanation generation').
prism_flag( error_on_cycle,oneof([on,off]), 'checks for cycles in derivations').
prism_flag( log_scale,  	oneof([on,off]), 'log scaling in inference and learning').
prism_flag( explicit_empty_expls, oneof([on,off]), 'affects structure of explanation graphs').
prism_flag( sort_hindsight,oneof([by_goal,by_prob]),   'how to sort subgoals in hindsight').

% affect viterbi algorithm
prism_flag( viterbi_mode,  oneof([ml,vb]),  'which parameters to use in Viterbi (probs or counts)').
prism_flag( log_viterbi,	oneof([on,off]), 'use log probabilities in Viterbi algorithm').
prism_flag( rerank, 			natural,         'affects Viterbi algorithm with hyperparameters').

% affect learning
prism_flag( learn_mode, 	oneof([ml, vb, both, ml_vt, vb_vt, both_vt]), 'switch between ML/EM and VB/EM'). 
prism_flag( data_source, 	oneof([data/1, file(_), none]), 'data source for learning').
prism_flag( init,     		oneof([none,random,noisy_u]), 'initialisation method for EM algorithm').
prism_flag( fix_init_order,oneof([on,off]), 'fixes switch initialisation order in EM algorithm').
prism_flag( std_ratio,  	nonneg,          'variance of noise when EM initialisation method is noisy_u').
prism_flag( restart,    	natural,         'number of restarts in EM algorithm').
prism_flag( epsilon,       nonneg,          'convergence threshold for EM algorirhm').
prism_flag( max_iterate,   natural,         'maximum iterations in EM algorithm (can be default)').
prism_flag( params_after_vbem, oneof([none,mean,max]),  'method of parameter estimation after VB/EM').
prism_flag( reset_hparams, oneof([on,off]), 'reset pseudocounts to default before VB/EM').

% DAEM learning
prism_flag( daem,     		oneof([on,off]), 'enable deterministic annealing EM algorithm').
prism_flag( itemp_init,    range(0,1),      'initial inverse temperature for DAEM algorithm').
prism_flag( itemp_rate,    range(1,inf),    'inverse temperature rate for DAEM algorithm').

prism_flag( mcmc_b,    	   natural,         'MCMC burn in').
prism_flag( mcmc_e,    	   natural,         'MCMC chain length').
prism_flag( mcmc_s,    	   natural,         'MCMC cycle length').

%% prism_flag_affects(-Name, -Subsystem) is nondet.
%
%  Relation between flag names and the subsytem affected, which
%  can be one of display, initialisation (of switches), explanation,
%  vitberbi, hindsight or learning.
prism_flag_affects( warn,	             display).
prism_flag_affects( verb,	             display).
prism_flag_affects( show_itemp,	       display).
prism_flag_affects( em_progress,        display).
prism_flag_affects( search_progress,	 display).
prism_flag_affects( learn_progress,	    display).
prism_flag_affects( mcmc_progress,	    display).
prism_flag_affects( write_call_events,  display).

prism_flag_affects( default_sw,           initialisation).
prism_flag_affects( default_sw_a,         initialisation).
prism_flag_affects( default_sw_d,         initialisation).

prism_flag_affects( clean_table,    explanation).
prism_flag_affects( reduce_copy,    explanation).
prism_flag_affects( log_scale,      explanation).
prism_flag_affects( error_on_cycle, explanation).
prism_flag_affects( sort_hindsight, hindsight).

prism_flag_affects( viterbi_mode,  viterbi).
prism_flag_affects( log_viterbi,   viterbi).
prism_flag_affects( rerank,        viterbi).
prism_flag_affects( log_scale,     viterbi).

prism_flag_affects( daem,           learning).
prism_flag_affects( itemp_init,     learning).
prism_flag_affects( itemp_rate,     learning).
prism_flag_affects( learn_mode, 	   learning).	
prism_flag_affects( data_source, 	learning).
prism_flag_affects( init,     		learning).
prism_flag_affects( fix_init_order, learning).
prism_flag_affects( std_ratio,     	learning).
prism_flag_affects( restart,    	   learning).
prism_flag_affects( epsilon,        learning).
prism_flag_affects( max_iterate,    learning).
prism_flag_affects( reset_hparams,  learning). 
prism_flag_affects( log_scale,      learning).
prism_flag_affects( mcmc_b,         learning).
prism_flag_affects( mcmc_e,         learning).
prism_flag_affects( mcmc_s,         learning).

% --------------------- INFORMATION ---------------------------

%% prism_show(+Subsystem:oneof([values,probs,counts,goals,flags])) is det.
%% prism_show(-Subsystem:oneof([values,probs,counts,goals,flags,stats])) is nondet.
%
%  Causes PRISM to print information to its output stream.
%  By default, this appears in a file 'prism.log' in the directory
%  where PRISM was started..
prism_show(values) :- prism_flush(show_values).
prism_show(probs)  :- prism_flush(show_sw).
prism_show(counts) :- prism_flush(show_sw_a).
prism_show(goals)  :- prism_flush(show_goals).
prism_show(flags)  :- prism_flush(show_prism_flags).
prism_show(stats)  :- prism_flush(prism_statistics).

prism_flush(G) :- prism(G), prism(flush_output).

%% prism_statistics(+StatName,-Value) is semidet.
%% prism_statistics(-StatName,-Value) is nondet.
%
%  Get values of various statistics on inference, learning, and
%  the explanation graph.
prism_statistics(S,V) :- prism_statistics(_,S,V).

%% prism_statistics(+Subsytem,+StatName,-Value) is semidet.
%% prism_statistics(-Subsystem,-StatName,-Value) is nondet.
%
%  Get values of various statistics on a PRISM subsytem, which can
%  be one of infer, learn or graph.
prism_statistics(infer,S,V) :- prism_nd(infer_statistics(S,V)).
prism_statistics(learn,S,V) :- prism_nd(learn_statistics(S,V)).
prism_statistics(graph,S,V) :- prism_nd(graph_statistics(S,V)).


% ------------------ SUPPORT PREDICATES -----------------------


interp(true).
interp(fail) :- !, fail.
interp(throw(Ex)) :- throw(Ex). %format('PRISM exception: ~w.\n',[Ex]).


send(Term) :-
	current_prism(_,Stream,_),
	write_canon(Stream,Term),
	write(Stream,'.\n'),
	flush_output(Stream).

recv(Term) :-
	current_prism(_,_,Stream),
	read_line_to_string(Stream, Input),
	parse(Term,Input).

parse(line(Input),Input).
parse(term(T,V),Input) :-
	(	catch(term_string(T,Input,[variable_names(V)]),_,fail) -> true
	;	(string_concat("query",_,Input) -> parse_bug(Input,T); true),
		format('PRISM> ~s\n',[Input]),fail).

parse_bug(Input,Term) :-
	term_string(Parse, Input, []),
	assert(parse_bug(Input,Parse,Term)).

compare(A,B) :- atomic(A), !, (A=B -> true; writeln(mismatch(A,B)), fail).
compare(A,B) :- 
	A =.. [FA|AA], B=.. [FB|BB], 
	compare(FA,FB),
	maplist(compare,AA,BB).
	
read_line(S,Line) :-
	wait_for_input([S],[_],0.001),
	read_line_to_codes(S,Line).

wait(TO) :-
	current_prism(_,_,Stream),
	wait_for_input([Stream],[_],TO).

% PRISM can't read 1e7. It wants 1.0e7.
% Therefore we have to reimplement write_canonical
write_canon(S,T) :- 
	(	numbervars(T,0,_,[singletons(true)]), 
		(canon(T,Codes,[]) -> true; throw(error(write_canon(T)))), 
		% format('~s',[Codes]), nl,
		format(S,'~s',[Codes]), fail
	;	true).

wrq(A,C1,C2) :- format(codes(C1,C2),'~q',[A]).

canon('$VAR'(N)) --> !, wrq('$VAR'(N)). 
canon([])        --> !, "[]".
canon([H|T])     --> !, "[", canon(H), comma_list(T), "]".
canon(X)         --> {float(X)}, !, {format(codes(C,T),'~15g',[X])}, float_codes(C,T).
canon(A)         --> {atomic(A)}, !, wrq(A).
canon(A)         --> {A=..[F,H|T]}, wrq(F), "(", canon(H), comma_list(T), ")".  

-(A,B,A,B).
float_codes(H,U)       --> {H==U}, !.
float_codes([46|T],U)  --> !, ".", T-U.
float_codes([101|T],U) --> !, ".0e", T-U.
float_codes([C|T],U)   --> [C], float_codes(T,U).

comma_list('$VAR'(N)) --> !, " | ", wrq('$VAR'(N)).
comma_list([])        --> !, {true}.
comma_list([H|T])     --> ", ", canon(H), comma_list(T). 

