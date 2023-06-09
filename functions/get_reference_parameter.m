function [reference_parameter] = get_reference_parameter(ID)
% returns reference parameters of given measurement ID
% Input:
%   ID: (1x1) ID of measurement as string
%
% Output:
%   struct reference_parameter with entries:
%   range:      (1xP) ranges of P people
%   doa:        (1xP) DoAs of P people
%   posture:    (1xP) Position of P people: 
%               sitting (0), laying (1), standing (2), sitting backwards to radar (3), 
%               standing sideways to the radar (4)
%   obstacle:   (1x1) Obstacle: free space (0), wall (1), door (2)
%   person_ID:  (1xP) ID of each person, which connects with the bitalino MAC IDs
%                   day 1: ["5A:E4", "5A:DF"], IDs: 1 - 11
%                   day 2: ["57:5F", "5A:E4", "57:60", "5B:3C", "5A:E3"], IDs: 12 - 37-2

ID_list = get_ID_list();

% ranges of persons from 1 to 3 meter
range = {[],[],[],...
         [1],[2],[3],[3],[1],[3.5],[2],[1.5],[2.5],[2.5],[1],... % 1 - 11
         [4,2,2,3,3], [2,2,2,2.75,2.5], [1.5,2,3,2.75,1], [2,1,3.5,2.5,2.5],... % 12 - 15, M13, M14 P4 sits closer to 2.75m than 3m
         [1.5,2,2,2,3], [1.5,1,1.5,2.5,2.5], [3,1.5,3,1,2.5], [1,3.5,2,2.5,1.5],... % 16 -19
         [1,2,2,2.5,3], [1,2,2.5,2.75,2], [1,3,3,2,1], [1,3.5,1.5,2.5,1.5],... % 20 - 23, M21 P4 sits closer to 2.75m than 3m
         [2.5,1,3,1.5,2], [2,2.5,2,1],[1,3,2,1.5],... % 24 - 26
         [3,4,3], [3,4,3], [3,4,3],... % 27
         [2.5,3,2.5], [2.5,3,2.5],[2.5,3,2.5],... % 28
         [2,1,2], [2,1,2], [2,1,2],... % 29 for M29 pdfs says 1.5m P3 sits always at 1m
         [2.5,2.5,2.5], [2.5,2.5,2.5], [2.5,2.5,2.5],... % 30
         [1.25,2.5,2], [1.25,2.5,2], [1.25,2.5,2],... % 31 P1 sits more at 1.25 in M31 and M31-1
         [3,2,2.5], [3,2,2.5], [3,2,2.5],... % 32
         [1.5,1.5], [1.5,1.5], [1.5,1.5],... % 33
         [2.5,2.5], [2.5,2.5], [2.5,2.5],... % 34
         [3,3], [3,3], [3,3],... % 35
         [3,3], [3,3], [3,3],... % 36
         [3,1.5], [3,1.5], [3,1.5]}.'; % 37

% angles of persons from -90 to 90 degree
doa = {[],[],[],...
         [-45],[-45],[0],[0],[0],[-45],[0],[-45],[0],[-45],[0],... % 1 - 11
         [0,-45,-60,-30,45], [-60,0,-45,30,45], [-45,30,-30,45,60], [-45,-30,0,45,60],... % 12 - 15
         [-60,0,-30,45,30], [-60,0,-30,30,45], [0,45,30,-60,-45], [-60,0,-45,30,45],... % 16 -19
         [-60,-30,-45,0,30], [-60,0,-30,30,45], [-60,30,0,45,60], [-60,0,-30,30,60],... % 20 - 23
         [-45,0,-30,30,60], [-30,30,0,60],[-60,30,0,60],... % 24 - 26
         [-30,0,30], [-30,0,30], [-30,0,30],... % 27
         [-45,30,-30], [-45,30,-30], [-45,30,-30]... % 28
         [-30,-60,0], [-30,-60,0], [-30,-60,0],... % 29
         [-45,0,-30], [-45,0,-30], [-45,0,-30],... % 30
         [60,45,-30], [60,45,-30], [60,45,-30],... % 31
         [0,30,-30], [0,30,-30], [0,30,-30],... % 32
         [-30,0], [-30,0], [-30,0],... % 33
         [-45,-30], [-45,-30], [-45,-30],... % 34
         [-30,0], [-30,0], [-30,0],... % 35
         [45,-30], [45,-30], [30,-30],... % 36, P1 sits in different place for 36 and 36-1
         [-60,-30], [-60,-30], [-45,-30]}.'; % 37, P1 sits in different place for 37-2

% sitting (0), laying (1), standing (2), sitting backwards to radar (3), standing sideways to the radar (4)
posture = {[],[],[],...
         [2],[1],[0],[0],[3],[0],[4],[0],[0],[2],[0],... % 1 - 11
         [0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,1], [0,0,1,0,0],... % 12 - 15
         [0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,0],... % 16 -19
         [0,0,0,1,0], [0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,1],... % 20 - 23
         [0,0,0,0,0], [0,0,0,1], [0,0,1,0],... % 24 - 26
         [0,0,0], [0,0,0], [0,0,0],... % 27
         [0,0,0], [0,0,0], [0,0,0],... % 28
         [0,0,0], [0,0,0], [0,0,0],... % 29
         [0,0,0], [0,0,0], [0,0,0],... % 30
         [0,0,0], [0,0,0], [0,0,0],... % 31
         [0,0,0], [0,0,0], [0,0,0],... % 32
         [0,0], [0,0], [0,0],... % 33
         [0,0], [0,0], [0,0],... % 34
         [0,0], [0,0], [0,0],... % 35
         [0,0], [0,0], [0,0],... % 36
         [0,0], [0,0], [0,0]}.'; % 37

% free space (0), wall (1), door (2)
obstacle = {[0],[1],[2],...
         [0],[0],[0],[2],[0],[0],[0],[2],[1],[1],[1],... % 1 - 11
         [0], [0], [0], [0],... % 12 - 15
         [0], [0], [0], [0],... % 16 -19
         [0], [0], [0], [0],... % 20 - 23
         [0], [0], [0],... % 24 - 26
         [0], [1], [2],... % 27
         [0], [1], [2],... % 28
         [0], [1], [2],... % 29
         [0], [1], [2],... % 30
         [0], [1], [2],... % 31
         [0], [1], [2],... % 32
         [0], [1], [2],... % 33
         [0], [1], [2],... % 34
         [0], [1], [2],... % 35
         [0], [1], [2],... % 36
         [0], [1], [2]}.'; % 37

% day 1: ["5A:E4", "5A:DF"], IDs: 1 - 11
% day 2: ["57:5F", "5A:E4", "57:60", "5B:3C", "5A:E3"], IDs: 12 - 37-2
% switched person 2 and 3
person_ID = {[],[],[],...
         [1],[1],[1],[1],[2],[2],[2],[2],[1],[1],[1],... % 1 - 11
         [1,2,3,4,5], [1,2,3,4,5], [1,2,3,4,5], [1,2,3,4,5],... % 12 - 15
         [1,2,3,4,5], [1,2,3,4,5], [1,2,3,4,5], [1,2,3,4,5],... % 16 -19
         [1,2,3,4,5], [1,2,3,4,5], [1,2,3,4,5], [1,2,3,4,5],... % 20 - 23
         [1,2,3,4,5], [1,2,3,4], [1,2,3,4],... % 24 - 26
         [1,2,3], [1,2,3], [1,2,3],... % 27
         [1,2,3], [1,2,3], [1,2,3],... % 28
         [1,2,3], [1,2,3], [1,2,3],... % 29
         [1,2,3], [1,2,3], [1,2,3],... % 30
         [1,2,3], [1,2,3], [1,2,3],... % 31
         [1,2,3], [1,2,3], [1,2,3],... % 32
         [1,3], [1,3], [1,3],... % 33
         [2,3], [2,3], [2,3],... % 34
         [2,3], [2,3], [2,3],... % 35
         [1,2], [1,2], [1,2],... % 36
         [1,2], [1,2], [1,2]}.'; % 37

ID_index = find(ID_list == ID);

if(isempty(ID_index))
    error("ID does not exist.")
elseif(length(ID_index) ~= 1)
    error("To many IDs, ID should be (1x1).")
end

reference_parameter.range = range{ID_index};
reference_parameter.doa = doa{ID_index};
reference_parameter.posture = posture{ID_index};
reference_parameter.obstacle = obstacle{ID_index};
reference_parameter.person_ID = person_ID{ID_index};



