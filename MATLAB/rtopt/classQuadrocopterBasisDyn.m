classdef(Abstract) classQuadrocopterBasisDyn <  classDyn & classTestEnv 
    %CLASSQUADROCOPTERBASISDYN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        F;
        J;
        H;
    end
    
    properties(Constant, GetAccess=private)
        %TODO: they come from the model -> model.getI, model.getI_M ...
        I = [0.0093886, 0.0093886, 0.018406];
        % Tr�gheitsmoment gesamt in kg m^2
        m   = 1.022;                                % Gesamtgewicht des Quadrokopters in kg
        I_M = 0.0001;     %ToDo                              % Tr�gheitsmoment der Motoren/Rotoren in kg m^2
        kT  = 1.5e-07;                              % N/(RPM^2)
        kQ  = 3e-09;                                % N*m/(RPM)^2
        d   = 0.22;                                 % Abstand in m
        g   = 9.81;
    end
        
    properties(GetAccess=private, SetAccess = protected)
         n_int;
         isEmptyF;
         isEmptyJ;
         isEmptyH;
    end
    
    methods(Access=private)
        function set_interval(obj, n_int)
            obj.n_int = n_int;
        end
    end
    
    methods
       
        function emptyResults(obj)
            
            %emptyResults@classModell(obj);
            
            obj.F   = [];
            obj.J  = [];
            obj.H = [];
            
            obj.isEmptyF = true;    
            obj.isEmptyJ = true;
            obj.isEmptyH = true;
        end

        function res = get.F(obj)
            if obj.isEmptyF 
                
                q   = obj.state(4:7    , :);
                v   = obj.state(8:10   , :);
                omega   = obj.state(11:13  , :);
                u   = obj.contr;
            
                Iges = obj.I;
                IM = obj.I_M;
                m = obj.m;
            
                kT  = obj.kT;
                kQ  = obj.kQ;
                d   = obj.d;
                g   = obj.g;
                
                t1 = q(4, :) .^ 2;
t2 = q(3, :) .^ 2;
t3 = q(1, :) .* q(4, :);
t4 = q(2, :) .* q(3, :);
t5 = q(1, :) .* q(3, :);
t6 = q(2, :) .* q(4, :);
t7 = q(2, :) .^ 2;
t8 = q(1, :) .* q(2, :);
t9 = q(3, :) .* q(4, :);
t10 = -2 .* t5 + 2 .* t6;
t11 = 2 .* t8 + 2 .* t9;
t12 = 1 - 2 .* t7 - 2 .* t2;
t13 = u(2, :) .^ 2;
t14 = u(1, :) .^ 2;
t15 = u(3, :) .^ 2;
t16 = u(4, :) .^ 2;
t17 = kT .* (t14 + t13 + t15 + t16) - m .* (omega(1, :) .* v(2, :) - omega(2, :) .* v(1, :) + t12 .* g);
t18 = -u(1, :) + u(2, :) - u(3, :) + u(4, :);
t19 = kT .* d;
t20 = omega(2, :) .* ((-Iges(3) + Iges(2)) .* omega(3, :) + IM .* t18) + t19 .* (t13 - t16);
t18 = omega(1, :) .* ((-Iges(1) + Iges(3)) .* omega(3, :) - IM .* t18) + t19 .* (t15 - t14);
cg = [(1 - 2 .* t2 - 2 .* t1) .* v(1, :) + (-2 .* t3 + 2 .* t4) .* v(2, :) + (2 .* t5 + 2 .* t6) .* v(3, :) (2 .* t3 + 2 .* t4) .* v(1, :) + (1 - 2 .* t7 - 2 .* t1) .* v(2, :) + (-2 .* t8 + 2 .* t9) .* v(3, :) t10 .* v(1, :) + t11 .* v(2, :) + t12 .* v(3, :) -(q(2, :) .* omega(1, :)) ./ 0.2e1 - (q(3, :) .* omega(2, :)) ./ 0.2e1 - (q(4, :) .* omega(3, :)) ./ 0.2e1 (q(1, :) .* omega(1, :)) ./ 0.2e1 + (q(3, :) .* omega(3, :)) ./ 0.2e1 - (q(4, :) .* omega(2, :)) ./ 0.2e1 (q(1, :) .* omega(2, :)) ./ 0.2e1 - (q(2, :) .* omega(3, :)) ./ 0.2e1 + (q(4, :) .* omega(1, :)) ./ 0.2e1 (q(1, :) .* omega(3, :)) ./ 0.2e1 + (q(2, :) .* omega(2, :)) ./ 0.2e1 - (q(3, :) .* omega(1, :)) ./ 0.2e1 -omega(2, :) .* v(3, :) + omega(3, :) .* v(2, :) - t10 .* g -omega(3, :) .* v(1, :) + omega(1, :) .* v(3, :) - t11 .* g 1 ./ m .* t17 1 ./ Iges(1) .* t20 1 ./ Iges(2) .* t18 1 ./ Iges(3) .* (omega(2, :) .* omega(1, :) .* (-Iges(2) + Iges(1)) + kQ .* (-t14 + t13 - t15 + t16))];


                obj.F = cg;
                obj.isEmptyF = false;
            end
            res = obj.F;
        end
        
        function res = get.J(obj)
            if obj.isEmptyJ
                q   = obj.state(4:7    , :);
                v   = obj.state(8:10   , :);
                omega   = obj.state(11:13  , :);
                u   = obj.contr;
            
                Iges = obj.I;
                IM = obj.I_M;
                m = obj.m;
            
                kT  = obj.kT;
                kQ  = obj.kQ;
                d   = obj.d;
                g   = obj.g;
                
                t228 = kT .* d;
t197 = q(3, :);
t227 = t197 .* g;
t198 = q(2, :);
t226 = t198 .* g;
t199 = q(1, :);
t225 = t199 .* g;
t224 = 1 - 2 .* t198 .^ 2;
t192 = u(4, :);
t193 = u(3, :);
t194 = u(2, :);
t195 = u(1, :);
t223 = (IM .* (t195 - t194 + t193 - t192));
t186 = Iges(3);
t180 = 1 ./ t186;
t222 = kQ .* t180;
t221 = -2 .* t228;
t220 = 2 .* t228;
t187 = Iges(2);
t181 = 0.1e1 ./ t187;
t185 = omega(1, :);
t219 = (t181 .* t185);
t188 = Iges(1);
t182 = 0.1e1 ./ t188;
t184 = omega(2, :);
t218 = (t182 .* t184);
t190 = v(2, :);
t196 = q(4, :);
t217 = t196 .* t190;
t191 = v(1, :);
t216 = t197 .* t191;
t189 = v(3, :);
t215 = t198 .* t189;
t214 = t199 .* t189;
t213 = t199 .* t190;
t212 = t199 .* t191;
t211 = t199 .* t196;
t210 = t199 .* t197;
t209 = t199 .* t198;
t208 = (-t186 + t187);
t207 = (-t188 + t186);
t206 = -2 .* t222;
t205 = 2 .* t222;
t204 = 2 ./ m .* kT;
t203 = IM .* t219;
t202 = IM .* t218;
t201 = (t180 .* (-t187 + t188));
t183 = omega(3, :);
t179 = t199 ./ 0.2e1;
t178 = -t198 ./ 0.2e1;
t177 = -t197 ./ 0.2e1;
t176 = -t196 ./ 0.2e1;
t174 = -2 .* t197 .^ 2;
t173 = -2 .* t196 .^ 2;
t172 = -t185 ./ 0.2e1;
t171 = t185 ./ 0.2e1;
t170 = -t184 ./ 0.2e1;
t169 = t184 ./ 0.2e1;
t168 = -t183 ./ 0.2e1;
t167 = t183 ./ 0.2e1;
t166 = -2 .* t226;
t165 = -2 .* t196 .* g;
t164 = IM .* t185;
t163 = IM .* t184;
t162 = t198 .* t197;
t161 = t198 .* t196;
t160 = t198 .* t191;
t159 = t198 .* t190;
t158 = t197 .* t196;
t157 = t197 .* t190;
t156 = t197 .* t189;
t155 = t196 .* t191;
t154 = t196 .* t189;
t1 = [1 4 -2 .* t217 + 2 .* t156; 1 5 2 .* t157 + 2 .* t154; 1 6 -4 .* t216 + 2 .* t159 + 2 .* t214; 1 7 -4 .* t155 - 2 .* t213 + 2 .* t215; 1 8 1 + t174 + t173; 1 9 -2 .* t211 + 2 .* t162; 1 10 2 .* t210 + 2 .* t161; 2 4 2 .* t155 - 2 .* t215; 2 5 2 .* t216 - 4 .* t159 - 2 .* t214; 2 6 2 .* t160 + 2 .* t154; 2 7 2 .* t212 - 4 .* t217 + 2 .* t156; 2 8 2 .* t211 + 2 .* t162; 2 9 t173 + t224; 2 10 -2 .* t209 + 2 .* t158; 3 4 -2 .* t216 + 2 .* t159; 3 5 2 .* t155 + 2 .* t213 - 4 .* t215; 3 6 -2 .* t212 + 2 .* t217 - 4 .* t156; 3 7 2 .* t160 + 2 .* t157; 3 8 -2 .* t210 + 2 .* t161; 3 9 2 .* t209 + 2 .* t158; 3 10 t174 + t224; 4 5 t172; 4 6 t170; 4 7 t168; 4 11 t178; 4 12 t177; 4 13 t176; 5 4 t171; 5 6 t167; 5 7 t170; 5 11 t179; 5 12 t176; 5 13 t197 ./ 0.2e1; 6 4 t169; 6 5 t168; 6 7 t171; 6 11 t196 ./ 0.2e1; 6 12 t179; 6 13 t178; 7 4 t167; 7 5 t169; 7 6 t172; 7 11 t177; 7 12 t198 ./ 0.2e1; 7 13 t179; 8 4 2 .* t227; 8 5 t165; 8 6 2 .* t225; 8 7 t166; 8 9 t183; 8 10 -t184; 8 12 -t189; 8 13 t190; 9 4 t166; 9 5 -2 .* t225; 9 6 t165; 9 7 -2 .* t227; 9 8 -t183; 9 10 t185; 9 11 t189; 9 13 -t191; 10 5 4 .* t226; 10 6 4 .* t227; 10 8 t184; 10 9 -t185; 10 11 -t190; 10 12 t191; 10 14 t195 .* t204; 10 15 t194 .* t204; 10 16 t193 .* t204; 10 17 t192 .* t204; 11 12 t182 .* (-t223 + t208 .* t183); 11 13 t208 .* t218; 11 14 -t202; 11 15 t182 .* ((t194 .* t220) + t163); 11 16 -t202; 11 17 t182 .* ((t192 .* t221) + t163); 12 11 t181 .* (t223 + t207 .* t183); 12 13 t207 .* t219; 12 14 t181 .* ((t195 .* t221) + t164); 12 15 -t203; 12 16 t181 .* ((t193 .* t220) + t164); 12 17 -t203; 13 11 t184 .* t201; 13 12 t185 .* t201; 13 14 t195 .* t206; 13 15 t194 .* t205; 13 16 t193 .* t206; 13 17 t192 .* t205;];

                
                obj.J = t1;
                obj.isEmptyJ = false;
            end
            res = obj.J;
        end
        
        function res = get.H(obj)
            if obj.isEmptyH
                q   = obj.state(4:7    , :);
                v   = obj.state(8:10   , :);
                omega   = obj.state(11:13  , :);
                u   = obj.contr;
                
             
                Iges = obj.I;
            
                IM = obj.I_M;
                m = obj.m;
            
                kT  = obj.kT;
                kQ  = obj.kQ;
                d   = obj.d;
                g   = obj.g;
                
                onesV = ones(1, obj.n_int); %Korregiere die Matrix Dimension

                t414 = g .* onesV;
t413 = kT .* onesV;
t389 = Iges(3);
t412 = onesV ./ t389;
t390 = Iges(2);
t385 = 1 ./ t390;
t411 = onesV .* t385;
t391 = Iges(1);
t386 = 1 ./ t391;
t410 = onesV .* t386;
t409 = v(3, :) .* onesV;
t408 = v(2, :) .* onesV;
t407 = v(1, :) .* onesV;
t406 = q(4, :) .* onesV;
t405 = q(3, :) .* onesV;
t404 = q(2, :) .* onesV;
t403 = q(1, :) .* onesV;
t402 = d .* t413;
t401 = kQ .* t412;
t400 = t385 .* t402;
t399 = t386 .* t402;
t388 = -onesV ./ 0.2e1;
t387 = onesV ./ 0.2e1;
t383 = -2 .* t414;
t382 = 2 .* t414;
t381 = 4 .* t414;
t380 = 2 ./ m .* t413;
t379 = -2 .* t403;
t378 = 2 .* t403;
t377 = -4 .* t404;
t376 = -2 .* t404;
t375 = 2 .* t404;
t374 = -4 .* t405;
t373 = -2 .* t405;
t372 = 2 .* t405;
t371 = -4 .* t406;
t370 = -2 .* t406;
t369 = 2 .* t406;
t368 = -4 .* t407;
t367 = -2 .* t407;
t366 = 2 .* t407;
t365 = -4 .* t408;
t364 = -2 .* t408;
t363 = 2 .* t408;
t362 = -4 .* t409;
t361 = -2 .* t409;
t360 = 2 .* t409;
t359 = IM .* t410;
t358 = IM .* t411;
t357 = -2 .* t401;
t356 = 2 .* t401;
t355 = (-t389 + t390) .* t410;
t354 = (-t391 + t389) .* t411;
t353 = (-t390 + t391) .* t412;
t1 = [1 4 6 t360; 1 4 7 t364; 1 4 9 t370; 1 4 10 t372; 1 5 6 t363; 1 5 7 t360; 1 5 9 t372; 1 5 10 t369; 1 6 4 t360; 1 6 5 t363; 1 6 6 t368; 1 6 8 t374; 1 6 9 t375; 1 6 10 t378; 1 7 4 t364; 1 7 5 t360; 1 7 7 t368; 1 7 8 t371; 1 7 9 t379; 1 7 10 t375; 1 8 6 t374; 1 8 7 t371; 1 9 4 t370; 1 9 5 t372; 1 9 6 t375; 1 9 7 t379; 1 10 4 t372; 1 10 5 t369; 1 10 6 t378; 1 10 7 t375; 2 4 5 t361; 2 4 7 t366; 2 4 8 t369; 2 4 10 t376; 2 5 4 t361; 2 5 5 t365; 2 5 6 t366; 2 5 8 t372; 2 5 9 t377; 2 5 10 t379; 2 6 5 t366; 2 6 7 t360; 2 6 8 t375; 2 6 10 t369; 2 7 4 t366; 2 7 6 t360; 2 7 7 t365; 2 7 8 t378; 2 7 9 t371; 2 7 10 t372; 2 8 4 t369; 2 8 5 t372; 2 8 6 t375; 2 8 7 t378; 2 9 5 t377; 2 9 7 t371; 2 10 4 t376; 2 10 5 t379; 2 10 6 t369; 2 10 7 t372; 3 4 5 t363; 3 4 6 t367; 3 4 8 t373; 3 4 9 t375; 3 5 4 t363; 3 5 5 t362; 3 5 7 t366; 3 5 8 t369; 3 5 9 t378; 3 5 10 t377; 3 6 4 t367; 3 6 6 t362; 3 6 7 t363; 3 6 8 t379; 3 6 9 t369; 3 6 10 t374; 3 7 5 t366; 3 7 6 t363; 3 7 8 t375; 3 7 9 t372; 3 8 4 t373; 3 8 5 t369; 3 8 6 t379; 3 8 7 t375; 3 9 4 t375; 3 9 5 t378; 3 9 6 t369; 3 9 7 t372; 3 10 5 t377; 3 10 6 t374; 4 5 11 t388; 4 6 12 t388; 4 7 13 t388; 4 11 5 t388; 4 12 6 t388; 4 13 7 t388; 5 4 11 t387; 5 6 13 t387; 5 7 12 t388; 5 11 4 t387; 5 12 7 t388; 5 13 6 t387; 6 4 12 t387; 6 5 13 t388; 6 7 11 t387; 6 11 7 t387; 6 12 4 t387; 6 13 5 t388; 7 4 13 t387; 7 5 12 t387; 7 6 11 t388; 7 11 6 t388; 7 12 5 t387; 7 13 4 t387; 8 4 6 t382; 8 5 7 t383; 8 6 4 t382; 8 7 5 t383; 8 9 13 onesV; 8 10 12 -onesV; 8 12 10 -onesV; 8 13 9 onesV; 9 4 5 t383; 9 5 4 t383; 9 6 7 t383; 9 7 6 t383; 9 8 13 -onesV; 9 10 11 onesV; 9 11 10 onesV; 9 13 8 -onesV; 10 5 5 t381; 10 6 6 t381; 10 8 12 onesV; 10 9 11 -onesV; 10 11 9 -onesV; 10 12 8 onesV; 10 14 14 t380; 10 15 15 t380; 10 16 16 t380; 10 17 17 t380; 11 12 13 t355; 11 12 14 -t359; 11 12 15 t359; 11 12 16 -t359; 11 12 17 t359; 11 13 12 t355; 11 14 12 -t359; 11 15 12 t359; 11 15 15 2 .* t399; 11 16 12 -t359; 11 17 12 t359; 11 17 17 -2 .* t399; 12 11 13 t354; 12 11 14 t358; 12 11 15 -t358; 12 11 16 t358; 12 11 17 -t358; 12 13 11 t354; 12 14 11 t358; 12 14 14 -2 .* t400; 12 15 11 -t358; 12 16 11 t358; 12 16 16 2 .* t400; 12 17 11 -t358; 13 11 12 t353; 13 12 11 t353; 13 14 14 t357; 13 15 15 t356; 13 16 16 t357; 13 17 17 t356;];

                
                obj.H = t1;
                obj.isEmptyH = false;
            end
            res = obj.H;
        end
    end
    
    methods(Test)
        function testFJH(obj)

            n_int_ = 3;
            n_state_ = obj.n_var;
            n_contr_ = obj.n_contr;
            
            obj.state = rand(n_state_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);

            F_ = obj.F;
            H_ = obj.H;
            J_ = obj.J;
        end
    end
    
end

