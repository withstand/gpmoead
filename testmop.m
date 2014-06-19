function mop = testmop( testname, dimension )
%Get test multi-objective problems from a given name.
%   The method get testing or benchmark problems for Multi-Objective
%   Optimization. The implemented problems included ZDT, OKA, KNO.
%   User get the corresponding test problem, which is an instance of class
%   mop, by passing the problem name and optional dimension parameters.

mop=struct('name',[],'nobj',[],'nvar',[],'domain',[],'func',[]);
switch lower(testname)
    case 'kno1'
        mop=kno1(mop);
    case 'sch2'
        mop = sch2(mop);
    case 'zdt1'
        mop=zdt1(mop, dimension);
    case 'kno3'
        mop = kno3(mop);
    case 'zdt2'
        mop = zdt2(mop,dimension);
    otherwise 
        error('Undefined test problem name');                
end 
end

%KNO1 function generator
function p=kno1(p)
 p.name='KNO1';
 p.nobj = 2;
 p.nvar = 2;
 p.domain= [0 0;3 3];
 p.func = vectorizefunction(@evaluate, p.nobj);
 
    %KNO1 evaluation function.
    function y = evaluate(x)
      y=zeros(2,1);
	  c = x(1)+x(2);
	  f = 9-(3*sin(2.5*c^0.5) + 3*sin(4*c) + 5 *sin(2*c+2));
	  g = (pi/2.0)*(x(1)-x(2)+3.0)/6.0;
	  y(1)= 20-(f*cos(g));
	  y(2)= 20-(f*sin(g)); 
    end
end

%KNO3 function generator
function p=kno3(p)
 p.name='KNO3';
 p.nobj = 3;
 p.nvar = 2;
 p.domain= [0 0;3 3];
 p.func = vectorizefunction(@evaluate, p.nobj);
 
    %KNO1 evaluation function.
    function y = evaluate(x)
      y=zeros(2,1);
	  c = x(1)+x(2);
	  f = 9-(3*sin(2.5*c^0.5) + 3*sin(4*c) + 5 *sin(2*c+2));
	  g = (pi/2.0)*(x(1)-x(2)+3.0)/6.0;
	  y(1)= 20-(f*cos(g));
	  y(2)= 20-(f*sin(g)); 
      y(3)= 20-(f*cos(g)*sin(g));
    end
end


%ZDT1 function generator
function p=zdt1(p,dim)
 p.name='ZDT1';
 p.nvar=dim;
 p.nobj=2;
 p.domain=[zeros(1,dim); ones(1,dim)];
 p.func=vectorizefunction(@evaluate, p.nobj);
 
    %KNO1 evaluation function.
    function y=evaluate(x)
        y=zeros(2,1);
        y(1) = x(1);
    	su = sum(x)-x(1);    
		g = 1 + 9 * su / (dim - 1);
		y(2) =g*(1 - sqrt(x(1) / g));
    end
end

function p=zdt2(p,dim)
 p.name='ZDT2';
 p.nvar=dim;
 p.nobj=2;
 p.domain=[zeros(1,dim); ones(1,dim)];
 p.func=vectorizefunction(@evaluate, p.nobj);
 
    %KNO1 evaluation function.
    function y=evaluate(x)
        y=zeros(2,1);
        y(1) = x(1);
    	su = sum(x(2:end));    
		g = 1 + 9 * su / (dim - 1);
		y(2) =g*(1 - (x(1) / g).^2);
    end
end

%ZDT1 function generator
function p=sch2(p)
 p.name='SCH2';
 p.nvar=1;
 p.nobj=2;
 p.domain=[-5 ;10];
 p.func=vectorizefunction(@evaluate, p.nobj);
 
    %KNO1 evaluation function.
    function y=evaluate(x)
        y=zeros(2,1);
        if x <= 1 
            y(1) = -x;
        elseif x <= 3
            y(1) = x-2;
        elseif x <=4
            y(1) = 4-x;
        else
            y(1) = x-4;
        end
        y(2) = (x-5)^2;
    end
end

 