function shiftArr = eoshift(a_Vector, a_Shift)

shiftArr = zeros(size(a_Vector));

if ( a_Shift < 0 )
    shiftArr(-a_Shift+1:length(shiftArr)) = a_Vector(1:length(a_Vector)+a_Shift);
elseif ( a_Shift > 0 )
    shiftArr(1:length(shiftArr)-a_Shift) = a_Vector(a_Shift+1:length(a_Vector));
elseif ( a_Shift == 0 )
    shiftArr = a_Vector;
end