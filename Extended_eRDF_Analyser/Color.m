function farbcode = Color(number)
% COLOR Liefert Farben und Mischungen davon als RGB-Tripel zwischen 0...1 zurück
% Eingegeben werden könnnen Zahlen von 1...13 wobei 0 schwarz und 14 weiß
% ist. Zwischen den einzelnen Farben wird per Linearkombination
% interpoliert, sodass z.b. die Werte 1.2, 1.5 und 1.8 Mischungen aus dem hellen und dunklen Blau ergeben.

if isnumeric(number)
    nummer=number;
else
    nummer=str2double(number);
end

if abs(nummer-round(nummer))<0.0001
    farbcode = ColorZuordnung(round(nummer));
else
    farbcode = ColorZuordnung(floor(nummer))*(1-nummer+floor(nummer)) + ColorZuordnung(ceil(nummer))*(1-ceil(nummer)+nummer);
end


end

function farbcode = ColorZuordnung(nummer)

% Liefert Farben als RGB-Tripel zwischen 0...1 zurück
switch nummer
    case 0
        farbe=[0,0,0];
    case 1
        farbe=[0,84,159];
    case 2
        farbe=[142,186,229];
    case 3
        farbe=[0,97,101];
    case 4
        farbe=[0,152,161];
    case 5
        farbe=[87,171,39];
    case 6
        farbe=[189,205,0];
    case 7
        farbe=[255,237,0];
    case 8
        farbe=[246,168,0];
    case 9
        farbe=[227,0,102];
    case 10
        farbe=[204,7,30];
    case 11
        farbe=[161,16,53];
    case 12
        farbe=[97,33,88];
    case 13
        farbe=[122,111,172];
    case 14
        farbe=[255,255,255];
    otherwise
        farbe=[0,0,0];
end
farbcode=farbe./255;
end