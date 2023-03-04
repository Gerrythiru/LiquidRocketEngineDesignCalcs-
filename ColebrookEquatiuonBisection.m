function [ColebrookFF] = ColebrookEquatiuonBisection(a, b, err, iter, Re, RelRough)
%This function calculates the colebrook-white friction factor by using the bisection method
%Inputs
%a Initial Lower Guess
%b Initial Higher Guess
%err error of the final result
%iter maximum number of iterations
%Re Reynolds number
%RelRough relative roughness of the tube (absolute roughness/inner diameter)
    
fa = -2*log10((RelRough/3.72)+(2.51/(Re*sqrt(a))))-(1/sqrt(a));   %Colebrook equation here
fb = -2*log10((RelRough/3.72)+(2.51/(Re*sqrt(b))))-(1/sqrt(b));   % in the form f(x) = 0

if fa > 0 || fb < 0
    error('variables not on opposite sides or correct order');
end

%fprintf('\n it \t\t a \t\t\t mid \t\t b \t\t F(mid)\t\t Error \n')

%loop to perform bisection
for i=1:iter
    
    %calculate mid and plug into function
    mid = (a+b)/2;
    fx = -2*log10((RelRough/3.72)+(2.51/(Re*sqrt(mid))))-(1/sqrt(mid));   %FUNCTION HERE
    MaxError = (b - a)/2;
    
    %compare to endpoints and replace
    if fx > 0
        b = mid;
    else 
        a = mid;
    end
    
    %fprintf('%.0f \t\t %.6f \t %.6f \t %.6f \t %.6f \t %.6f \n', i, a, mid, b, fx, MaxError)
    if MaxError < err
        %fprintf('Converges at %.6f\n',  mid);
        break;
    end
end
ColebrookFF = mid;
end