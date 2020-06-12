classdef intvclass2 <handle 
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
%%    
    properties
    end
    
    methods
%%        
        function obj = intvclass2()
        end
%% Interval Definition   
        function res = def(obj,varargin)
            count=(nargin-1);

            if iscell(varargin{1})
                varargin=varargin{1};
                count=count+1;               
            end    

            if count>2; error('Interval should be defined by only 2 values'); end
            
            if count==1
                n1=varargin{1};
                res = obj.def(n1,n1);
            else
                n1=varargin{1};
                n2=varargin{2};
                res = [{min(n1,n2)} {max(n1,n2)}];
            end
        end

        
%% Is Interval? check
        function res = isinterval(obj,x)
            res=0;
            if iscell(x)==1
                if length(x)==2
                    if length(x{1})==length(x{2})
                        if unique(le(x{1},x{2}))==1
                            res=1;
                        else
                            warning('Interval elements are not arranged correctly')
                        end
                    else
                        warning('Number of interval elements are not the same')                        
                    end
                end
            end
        end

%% Thin Interval check
%Returns 1 if the input is a single number, or an interval is thin
        function res = thin(obj,x)
            if iscell(x)==1
                if length(x)==2
                    if isequal(x{1},x{2})
                        res=1;
                    else
                        res=0;
                    end
                else
                    error('only one cell in interval')
                end
            else
                res=1;
            end
        end  
%% Interval to Matrix
        function res = int2mat(obj,x)
            if ~obj.isinterval(x); x=obj.def(x); end
            res= [x{1} x{2}];
        end
        
%% Union
        function res = union(obj,x,y) %needs revision
            res = [{min(x{1},y{1})} {max(x{2},y{2})}];
        end     
        
%% Intersection
        function res = intersec(obj,x,y) %needs revision
            res = [{max(x{1},y{1})} {min(x{2},y{2})}];
        end  
%% Addition       
        function res = add(obj,x,y)
            if ~obj.isinterval(x); x=obj.def(x); end
            if ~obj.isinterval(y); y=obj.def(y); end
            if length(x{1})==length(y{1})
                res{1} = x{1}+ y{1};
                res{2} = x{2}+ y{2};
            else
                error('Intervals size does not match')
            end
        end
%% Subtraction       
        function res = sub(obj,x,y)
            if ~obj.isinterval(x); x=obj.def(x); end
            if ~obj.isinterval(y); y=obj.def(y); end
            y=obj.def(-y{1},-y{2});
            res = obj.add(x,y);
        end
%% Element Multiplication        
        function res = elemult(obj,x,y)
            
            if 0<=x(1) && 0<=y(1)    
                res(1) = x(1)* y(1);
                res(2) = x(2)* y(2);
                return;
            end
            if x(1)<0 && 0<x(2) && 0<=y(1)    
                res(1) = x(1)* y(2);
                res(2) = x(2)* y(2);
                return;
            end
            if x(2)<=0 && 0<=y(1)    
                res(1) = x(1)* y(2);
                res(2) = x(2)* y(1);
                return;
            end
            if 0<=x(1) && y(1)<0 && y(2)>0    
                res(1) = x(2)* y(1);
                res(2) = x(2)* y(2);
                return;
            end
            if x(2)<=0 && y(1)<0 && y(2)>0    
                res(1) = x(1)* y(2);
                res(2) = x(1)* y(1);
                return;
            end 
            if 0<=x(1) && y(2)<=0    
                res(1) = x(2)* y(1);
                res(2) = x(1)* y(2);
                return;
            end
            if x(1)<0 && x(2)>0 && y(2)<=0    
                res(1) = x(2)* y(1);
                res(2) = x(1)* y(1);
                return;
            end  
            if x(2)<=0 && y(2)<=0    
                res(1) = x(2)* y(2);
                res(2) = x(1)* y(1);
                return;
            end  
            if x(1)<0 && x(2)>0 && y(1)<0 && y(2)>0
                temp1=x(1)*y(2);temp2=x(2)*y(1);
                temp3=x(1)*y(1);temp4=x(2)*y(2);
                res(1) = min(temp1,temp2);
                res(2) = max(temp3,temp4);
                return;
            end                      
        end
%% Multiplication
        function res = mult(obj,A,B)
            
            if ~obj.isinterval(A); A=obj.def(A); end
            if ~obj.isinterval(B); B=obj.def(B); end
            
            [rowsA, colA] = size(A{1});
            [rowsB, colB] = size(B{1});
            
            res=obj.def(zeros(rowsA,colB));
            
            for row=1:rowsA
                for col=1:colB
                    sum=obj.def(0);
                    for k=1:colA
                        x=[A{1}(row,k) A{2}(row,k)];
                        y=[B{1}(k,col) B{2}(k,col)];
                        prod=obj.elemult(x,y);
                        prod=obj.def(prod(1),prod(2));
                        sum=obj.add(sum,prod);
                    end
                    res{1}(row,col)=sum{1};
                    res{2}(row,col)=sum{2};
                end
            end
        end
        
%% Division
        function res = div(obj,x,y)
            if ~obj.isinterval(x); x=obj.def(x); end
            if ~obj.isinterval(y); y=obj.def(y); end
                if all(unique(sign(y{1})==sign(y{2}))) ...
                        && all(unique(y{1}~=0)) ...
                        && all(unique(y{2}~=0))
                    invy=obj.def(y{1}.^(-1),y{2}.^(-1));
                    res = obj.mult(x,invy);
                else
                    res=[];
                    error('Undefined operation: Division by interval containing zero. ')
                end
        end
%% Power        
%         function res = pow(obj,x,power)
% %           warning('power should be integer')
%             if obj.thin(x); res=x^power; return; end
%             if x==0; res=0; return; end
%             if power==0; res=obj.def(1,1); return; end
%             if power==1; res=x; return; end
%             res=x;
%             for i=2:power
%                 res=obj.mult(res,x);
%             end
%         end

                    
            

    end
end
