% bprepl : evaluator service to run in B-Prolog
%
% Sits on standard input reading terms and sending replies
% Query must be of the form:
%    query(+Id,+Query,+Vars)
% The service calls Query and replies with
%    query(+Id,+Vars,+Reply)
% Reply can be true, fail, or throw(Exception).
% Queries must be semideterministic as multiple solutions
% are not suppoerted at present.

:- dynamic logging/0.

repl(In,Out) :-
	once((repeat,read_valid_term(In,Term,Vars))),
	(logging, write('-> '), labelvars(Vars), reply_term(Term), fail; true),
	(	Term=end_of_file -> true
	;	interp(Term,Reply),
		labelvars(Vars),
		(logging -> write('<- '), reply_term(Reply); true),
		reply_term(Out,Reply),
		repl(In,Out)
	).

read_valid_term(Stream,Term,Vars) :-
	catch( read_term(Stream,Term,[variable_names(Vars)]),
		Ex,(writeln(Ex),fail)).

reply_term(T) :- current_output(S), reply_term(S,T).
reply_term(S,T) :- write_canon(S,T), write(S,'.\n'). 

% interp is not allowed to fail or throw and exception
interp(query(Id,Query,Vars),query(Id,Vars,Reply)) :- !,
	catch((call(Query) -> Reply=true; Reply=fail), E, Reply=throw(E)).
interp(_,unrecognised).

% this is so variables in reply match names of variables in query
labelvars([N=V|T]) :- (var(V) -> V='$VAR'(N); true), labelvars(T).
labelvars([]).

% built in write_canonical doesn't handle $VAR terms well
write_canon(S,'$VAR'(N)) :- !, write(S,N).
write_canon(S,[]) :-!, write(S,[]).
write_canon(S,[H|T]) :- !, 
	write(S,'['), write_canon(S,H),
	write_list_tail(S,T),
	write(S,']').

write_canon(S,A) :- float(A), !, format(S,'~15g',[A]).
write_canon(S,A) :- atomic(A), !, writeq(S,A).
write_canon(S,A) :- A=..[F,H|T], writeq(S,F), 
	write(S,'('), write_canon(S,H),
	write_list_tail(S,T),
	write(S,')').

write_list_tail(_,[]).
write_list_tail(S,[H|T]) :- 
	write(S,', '), write_canon(S,H), 
	write_list_tail(S,T).


open_log(file(File),Con) :- open(File,append,Con).

main :- main(file('bprepl.log')).

main(LogSpec) :- 
	open_log(LogSpec,Con),
	set_output(Con),
	writeln('-- starting'),
	catch((
		reply_term(user_output,repl(ready)),
		repl(user_input,user_output)),
		Ex,(write('** exception: '), writeln(Ex))),
	writeln('-- terminating'),
	set_output(user_output).

log(on) :- logging -> true; assert(logging).
log(off) :- retractall(logging).

% PRISM stuff

:- dynamic dynamic_values/2.

set_dynamic_values(I,V) :-
	(	dynamic_values(I,V2)
	->	(	V=V2 -> true
		; 	retractall(dynamic_values(I,_)), 
			assert(dynamic_values(I,V))
		)
	;	assert(dynamic_values(I,V))
	).


