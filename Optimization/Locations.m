function [x,y] = Locations(Location_Flag, n_wec)

msg = 'Array is not prescribed for this number of WECs.';

switch upper(Location_Flag)

    case 'NEAR-SYMMETRICAL'

        if n_wec == 1

            x = 0; y = 0;

        elseif n_wec == 3

            % WEC x location
            x                  = [0; 0; 0; 0; 0];

            % WEC y location
            y                  = [0; 0; 0; 0; 0];

        elseif n_wec == 5
            % WEC x location
            x                  = [0; 30*cos(pi/6); 30*cos(pi/6); 60*cos(pi/9); 60*cos(pi/9)];

            % WEC y location
            y                  = [0; 30*sin(pi/6); -30*sin(pi/6); 60*sin(pi/9); -60*sin(pi/9)];
        elseif n_wec == 6
            % WEC x location
            x                  = [0; 30*cos(pi/6); 30*cos(pi/6); 60*cos(pi/9); 60*cos(pi/9); 80];

            % WEC y location
            y                  = [0; 30*sin(pi/6); -30*sin(pi/6); 60*sin(pi/9); -60*sin(pi/9); 0];

        else
            error(msg)
        end

    case 'FAR-SYMMETRICAL'

        if n_wec == 1
            
            x = 0; y = 0;

        elseif n_wec == 5
            % WEC x location
            x                  = [0; 600*cos(pi/3.8); 600*cos(pi/3.8); 1200*cos(pi/3.8); 1200*cos(pi/3.8)];

            % WEC y location
            y                  = [0; 600*sin(pi/3.8); -600*sin(pi/3.8); 1200*sin(pi/3.8); -1200*sin(pi/3.8)];
        else
            error(msg)
        end

    case 'ASYMMETRICAL'

        if n_wec == 1

            x = 0; y = 0;

        elseif n_wec == 5
            % WEC x location
            x                  = [0; 30*cos(pi/3.14); 25*cos(pi/6); 35*cos(pi/6); 70*cos(pi/2.2)];

            % WEC y location
            y                  = [0; 30*sin(pi/6); -45*sin(pi/3); 30*sin(pi/1.5); -65*sin(pi/2)];
        else
            error(msg)
        end

    case 'COLUMN'

        if n_wec == 1

            x = 0; y = 0;

        elseif n_wec == 5
            % WEC x location
            x                  = [0; 30*cos(pi/2); 25*cos(pi/2); 35*cos(pi/2); 70*cos(pi/2)];

            % WEC y location
            y                  = [0; 30*sin(pi/2); -30*sin(pi/2); 60*sin(pi/2); -60*sin(pi/2)];
        else
            error(msg)
        end

    case 'ROW'

        if n_wec == 1

            x = 0; y = 0;

        elseif n_wec == 5
            % WEC x location
            x                  = [0; 30*sin(pi/2); 60*sin(pi/2); 90*sin(pi/2); 120*sin(pi/2)];


            % WEC y location
            y                  = [0; 30*cos(pi/2); 25*cos(pi/2); 35*cos(pi/2); 70*cos(pi/2)];
        else
            error(msg)
        end

end
end
