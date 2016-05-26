function [S, Sed] = getStoichMatrix()
	S=[-1 1 0 0 0 0;1 -1 0 0 0 0;0 0 1 -1 0 0;0 0 0 0 1 -1];
	Sed=[-1 0 -1 0 0 0;0 -1 0 0 0 0;0 0 0 -1 -1 0;0 0 0 0 0 -1];
end