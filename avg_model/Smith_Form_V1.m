


%A=input('Enter your matrix(n*n):>')
Ehsan_Matrix_Dimension=size(A);

while (Ehsan_Matrix_Dimension(1,1) ~= Ehsan_Matrix_Dimension(1,2))
    clear A;
    A=input('Size Mismatch on matrix dimension.Enter your matrix again(n*n):>');
end

Ehsan_Matrix_Dimension=Ehsan_Matrix_Dimension(1,1);

Ehsan_Temp_Level_Counter_i=0;
Ehsan_Change_Available=0;
Ehsan_Temporary_Indicator=0;
Ehsan_Level_Of_Operation=1;
Ehsan_Operation_Failar = 0;
Ehsan_Phase_One_Finished = 0;
Ehsan_Operation_Finished = 0;

syms t;

for Ehsan_i_Temporary=1:Ehsan_Matrix_Dimension
    for Ehsan_j_Temporary=1:Ehsan_Matrix_Dimension
        if Ehsan_i_Temporary==Ehsan_j_Temporary
            I(Ehsan_i_Temporary,Ehsan_j_Temporary)=1;
        else
            I(Ehsan_i_Temporary,Ehsan_j_Temporary)=0;
        end
    end
end
        
B=t.*I-A;
%Bikhodi ast engar.Ehsan_Level_Of_Operation=Ehsan_Temp_Level_Of_Operation;

while Ehsan_Level_Of_Operation < Ehsan_Matrix_Dimension

'Level_Of_Operation='
Ehsan_Level_Of_Operation

Ehsan_Sum_Of_Elements_In_First_Coloumn = 0;
Ehsan_Sum_Of_Elements_In_First_Row = 0;



%
for Ehsan2_i_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Matrix_Dimension
    for Ehsan2_j_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Matrix_Dimension
        Ehsan2_Temporary_Ranks_Of_Elements = sym2poly(expand(B(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter)));
        Ehsan2_Polynomials_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary_Ranks_Of_Elements(1,1);
        Ehsan2_Temporary2_Ranks_Of_Elements = size(Ehsan2_Temporary_Ranks_Of_Elements);
        Ehsan2_Ranks_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary2_Ranks_Of_Elements(1,2);
    end
end

clear Ehsan2_i_Temporary_Counter;
clear Ehsan2_Temporary_Ranks_Of_Elements;
clear Ehsan2_Temporary2_Ranks_Of_Elements;
Ehsan_Sum_Of_Elements_In_First_Coloumn = 0;
for Ehsan_Temporary_Counter_n = Ehsan_Level_Of_Operation + 1:Ehsan_Matrix_Dimension
    Ehsan_Sum_Of_Elements_In_First_Coloumn=Ehsan_Sum_Of_Elements_In_First_Coloumn+abs(B(Ehsan_Temporary_Counter_n,Ehsan_Level_Of_Operation));
end
clear Ehsan_Temporary_Counter_n;

if Ehsan_Sum_Of_Elements_In_First_Coloumn ~= 0
    Ehsan_All_Are_Resorted = 0;
    Ehsan_Temporary_Counter_i = 1;
    Ehsan_Temporary_element = Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation + Ehsan_Temporary_Counter_i,Ehsan_Level_Of_Operation);
    
    %chek mikonim ke avalin elemani ke migozarim gheyre sefr bashad.%
    while (expand(B(Ehsan_Level_Of_Operation + Ehsan_Temporary_Counter_i,Ehsan_Level_Of_Operation))) == 0
        Ehsan_Temporary_Counter_i = Ehsan_Temporary_Counter_i + 1;
        Ehsan_Temporary_element = Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation + Ehsan_Temporary_Counter_i,Ehsan_Level_Of_Operation);
    end
    
    Ehsan_Lowest_i = Ehsan_Temporary_Counter_i;
    
%bikhodi    if (Ehsan_Level_Of_Operation +1 + Ehsan_Temporary_Counter_i) <= Ehsan_Matrix_Dimension
        while (Ehsan_All_Are_Resorted == 0) && ((Ehsan_Level_Of_Operation +1 + Ehsan_Temporary_Counter_i) <= Ehsan_Matrix_Dimension)
            if  Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation +1 + Ehsan_Temporary_Counter_i,Ehsan_Level_Of_Operation) < Ehsan_Temporary_element && (expand(B(Ehsan_Level_Of_Operation + Ehsan_Temporary_Counter_i + 1,Ehsan_Level_Of_Operation))) ~= 0
                Ehsan_Temporary_element = Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation +1 + Ehsan_Temporary_Counter_i,Ehsan_Level_Of_Operation);
                Ehsan_Lowest_i = Ehsan_Temporary_Counter_i;
            end
        Ehsan_Temporary_Counter_i = Ehsan_Temporary_Counter_i + 1;
        
        if Ehsan_Temporary_Counter_i == (Ehsan_Matrix_Dimension - Ehsan_Level_Of_Operation)
            Ehsan_All_Are_Resorted=1;
            
        end
        end
%bikhodi    end
    if Ehsan_Temporary_Counter_i == (Ehsan_Matrix_Dimension - Ehsan_Level_Of_Operation)
         Ehsan_All_Are_Resorted=1;
            
    end
    if Ehsan_All_Are_Resorted == 1
        C=B;expand(C);
        C(Ehsan_Level_Of_Operation + 1,:)=B(Ehsan_Level_Of_Operation + Ehsan_Lowest_i,:);
        C(Ehsan_Level_Of_Operation + Ehsan_Lowest_i,:)=B(Ehsan_Level_Of_Operation + 1,:);
        B=C;expand(B);
    
        for Ehsan2_i_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Matrix_Dimension
            for Ehsan2_j_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Matrix_Dimension
                Ehsan2_Temporary_Ranks_Of_Elements = sym2poly(expand(B(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter)));
                Ehsan2_Polynomials_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary_Ranks_Of_Elements(1,1);
                Ehsan2_Temporary2_Ranks_Of_Elements = size(Ehsan2_Temporary_Ranks_Of_Elements);
                Ehsan2_Ranks_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary2_Ranks_Of_Elements(1,2);
            end
        end
    end
    
clear C;
clear Ehsan_Temporary_Counter_i;
clear Ehsan_Lowest_i;
clear Ehsan_Temporary_element;
clear Ehsan_All_Are_Resorted;
end



Ehsan_Operation_Finished=0;

if Ehsan_Sum_Of_Elements_In_First_Coloumn ~=0
    while Ehsan_Operation_Finished==0 && Ehsan_Operation_Failar == 0
    
        if (Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation) <= Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation + 1,Ehsan_Level_Of_Operation)) && isreal(expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation))) == 0 && B(Ehsan_Level_Of_Operation+1,Ehsan_Level_Of_Operation) ~= 0
            B(Ehsan_Level_Of_Operation+1,:)=expand(expand(B(Ehsan_Level_Of_Operation,:)).*(((-Ehsan2_Polynomials_Of_Elements(Ehsan_Level_Of_Operation + 1,Ehsan_Level_Of_Operation)/Ehsan2_Polynomials_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation)))*t^((Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation + 1,Ehsan_Level_Of_Operation)-Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation)))))+expand((B(Ehsan_Level_Of_Operation+1,:)));expand(B);
        
            for Ehsan2_i_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Level_Of_Operation+1
                for Ehsan2_j_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Matrix_Dimension
                    Ehsan2_Temporary_Ranks_Of_Elements = sym2poly(expand(B(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter)));
                    Ehsan2_Polynomials_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary_Ranks_Of_Elements(1,1);
                    Ehsan2_Temporary2_Ranks_Of_Elements = size(Ehsan2_Temporary_Ranks_Of_Elements);
                    Ehsan2_Ranks_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary2_Ranks_Of_Elements(1,2);
                end
            end

        clear Ehsan2_i_Temporary_Counter;
        clear Ehsan2_Temporary_Ranks_Of_Elements;
        clear Ehsan2_Temporary2_Ranks_Of_Elements;
        end

        if (Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation) > Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation + 1,Ehsan_Level_Of_Operation)) && isreal(expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation))) == 0 && B(Ehsan_Level_Of_Operation+1,Ehsan_Level_Of_Operation) ~= 0
            B(Ehsan_Level_Of_Operation,:)=((expand(B(Ehsan_Level_Of_Operation+1,:))).*(expand(((-Ehsan2_Polynomials_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation)/Ehsan2_Polynomials_Of_Elements(Ehsan_Level_Of_Operation+1,Ehsan_Level_Of_Operation)))*(t^(((Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation)-Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation+1,Ehsan_Level_Of_Operation))))))))+expand(B(Ehsan_Level_Of_Operation,:));expand(B);
            
            for Ehsan2_i_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Level_Of_Operation+1
                for Ehsan2_j_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Matrix_Dimension
                    Ehsan2_Temporary_Ranks_Of_Elements = sym2poly(expand(B(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter)));
                    Ehsan2_Polynomials_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary_Ranks_Of_Elements(1,1);
                    Ehsan2_Temporary2_Ranks_Of_Elements = size(Ehsan2_Temporary_Ranks_Of_Elements);
                    Ehsan2_Ranks_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary2_Ranks_Of_Elements(1,2);
                end
            end

            clear Ehsan2_i_Temporary_Counter;
            clear Ehsan2_Temporary_Ranks_Of_Elements;
            clear Ehsan2_Temporary2_Ranks_Of_Elements;
            if isreal(expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation))) == 1
                Ehsan_Operation_Finished=1;
            end
        end
        if isreal(expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation))) == 1
            Ehsan_Operation_Finished=1;
        end
        if B(Ehsan_Level_Of_Operation + 1,Ehsan_Level_Of_Operation) == 0
            Ehsan_Operation_Failar = 1;
        end

    end

end




if Ehsan_Sum_Of_Elements_In_First_Coloumn == 0
    Ehsan_Sum_Of_Elements_In_First_Row = 0 ;
    for Ehsan_Temporary_Counter_n2 = Ehsan_Level_Of_Operation + 1:Ehsan_Matrix_Dimension
        Ehsan_Sum_Of_Elements_In_First_Row=Ehsan_Sum_Of_Elements_In_First_Row+abs(B(Ehsan_Level_Of_Operation,Ehsan_Temporary_Counter_n2));
    end
    clear Ehsan_Temporary_Counter_n2;
    
    if Ehsan_Sum_Of_Elements_In_First_Row ~=0
        
        %ghesmate aval ke avaz kardane jaye sotoonha mishavad dar inja ra
        %naneveshtei,mesle ghabl amal kon agar L+1 sefr bood ba yek gheyre
        %sefr bayad avaz shavad.line 65
        Ehsan_All_Are_Resorted = 0;
        Ehsan_Temporary_Counter_i = 1;
        Ehsan_Temporary_element = Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation + Ehsan_Temporary_Counter_i,Ehsan_Level_Of_Operation);
    
        %chek mikonim ke avalin elemani ke migozarim gheyre sefr bashad.%
        while (expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + Ehsan_Temporary_Counter_i))) == 0
            Ehsan_Temporary_Counter_i = Ehsan_Temporary_Counter_i + 1;
            Ehsan_Temporary_element = Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + Ehsan_Temporary_Counter_i);
        end
    
        Ehsan_Lowest_i = Ehsan_Temporary_Counter_i;
    
      %bikhodi    if (Ehsan_Level_Of_Operation +1 + Ehsan_Temporary_Counter_i) <= Ehsan_Matrix_Dimension
        while (Ehsan_All_Are_Resorted == 0) && ((Ehsan_Level_Of_Operation +1 + Ehsan_Temporary_Counter_i) <= Ehsan_Matrix_Dimension)
            if  Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation +1 + Ehsan_Temporary_Counter_i) < Ehsan_Temporary_element && (expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + Ehsan_Temporary_Counter_i + 1))) ~= 0
                Ehsan_Temporary_element = Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation +1 + Ehsan_Temporary_Counter_i);
                Ehsan_Lowest_i = Ehsan_Temporary_Counter_i;
            end
            Ehsan_Temporary_Counter_i = Ehsan_Temporary_Counter_i + 1;
        
            if Ehsan_Temporary_Counter_i == (Ehsan_Matrix_Dimension - Ehsan_Level_Of_Operation)
                Ehsan_All_Are_Resorted=1;
            
            end
        end
      %bikhodi    end
    if Ehsan_Temporary_Counter_i == (Ehsan_Matrix_Dimension - Ehsan_Level_Of_Operation)
        Ehsan_All_Are_Resorted=1;
    end
    
        if Ehsan_All_Are_Resorted == 1
            C=B;expand(C);
            C(:,Ehsan_Level_Of_Operation + 1)=B(:,Ehsan_Level_Of_Operation + Ehsan_Lowest_i);
            C(:,Ehsan_Level_Of_Operation + Ehsan_Lowest_i)=B(:,Ehsan_Level_Of_Operation + 1);
            B=C;expand(B);
    
            for Ehsan2_i_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Matrix_Dimension
                for Ehsan2_j_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Matrix_Dimension
                    Ehsan2_Temporary_Ranks_Of_Elements = sym2poly(expand(B(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter)));
                    Ehsan2_Polynomials_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary_Ranks_Of_Elements(1,1);
                    Ehsan2_Temporary2_Ranks_Of_Elements = size(Ehsan2_Temporary_Ranks_Of_Elements);
                    Ehsan2_Ranks_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary2_Ranks_Of_Elements(1,2);
                end
            end
        end
    
    clear C;
    clear Ehsan_Temporary_Counter_i;
    clear Ehsan_Lowest_i;
    clear Ehsan_Temporary_element;
    clear Ehsan_All_Are_Resorted;

        %amaliat mesle sotooni%
        while Ehsan_Operation_Finished==0 && Ehsan_Operation_Failar == 0
            if (Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation) <= Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + 1)) && isreal(expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation))) == 0 && B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + 1) ~= 0
                B(:,Ehsan_Level_Of_Operation+1)=expand(expand(B(:,Ehsan_Level_Of_Operation)).*(((-Ehsan2_Polynomials_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + 1)/Ehsan2_Polynomials_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation)))*t^((Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + 1)-Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation)))))+expand(B(:,Ehsan_Level_Of_Operation+1));expand(B);
        
                for Ehsan2_i_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Matrix_Dimension
                    for Ehsan2_j_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Level_Of_Operation + 1
                        Ehsan2_Temporary_Ranks_Of_Elements = sym2poly(expand(B(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter)));
                        Ehsan2_Polynomials_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary_Ranks_Of_Elements(1,1);
                        Ehsan2_Temporary2_Ranks_Of_Elements = size(Ehsan2_Temporary_Ranks_Of_Elements);
                        Ehsan2_Ranks_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary2_Ranks_Of_Elements(1,2);
                    end
                end

            clear Ehsan2_i_Temporary_Counter;
            clear Ehsan2_Temporary_Ranks_Of_Elements;
            clear Ehsan2_Temporary2_Ranks_Of_Elements;
            end

            if (Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation) > Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + 1)) && isreal(expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation))) == 0 && B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + 1) ~= 0
                B(:,Ehsan_Level_Of_Operation)=expand((expand(B(:,Ehsan_Level_Of_Operation+1))).*((-(Ehsan2_Polynomials_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation)/(Ehsan2_Polynomials_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + 1))))*(t^(((Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation)-Ehsan2_Ranks_Of_Elements(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + 1))))))+B(:,Ehsan_Level_Of_Operation));expand(B);
            
                for Ehsan2_i_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Matrix_Dimension
                    for Ehsan2_j_Temporary_Counter = Ehsan_Level_Of_Operation:Ehsan_Level_Of_Operation+1
                        Ehsan2_Temporary_Ranks_Of_Elements = sym2poly(expand(B(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter)));
                        Ehsan2_Polynomials_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary_Ranks_Of_Elements(1,1);
                        Ehsan2_Temporary2_Ranks_Of_Elements = size(Ehsan2_Temporary_Ranks_Of_Elements);
                        Ehsan2_Ranks_Of_Elements(Ehsan2_i_Temporary_Counter,Ehsan2_j_Temporary_Counter) = Ehsan2_Temporary2_Ranks_Of_Elements(1,2);
                    end
                end

                clear Ehsan2_i_Temporary_Counter;
                clear Ehsan2_Temporary_Ranks_Of_Elements;
                clear Ehsan2_Temporary2_Ranks_Of_Elements;
                if isreal(expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation))) == 1
                    Ehsan_Operation_Finished=1;
                end
            end
            if expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation + 1)) == 0 && isreal(expand(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation))) == 0
                Ehsan_Operation_Failar = 1;
            end

        end

    else if Ehsan_Sum_Of_Elements_In_First_Coloumn == 0 && Ehsan_Sum_Of_Elements_In_First_Row == 0 && Ehsan_Level_Of_Operation < Ehsan_Matrix_Dimension
            %levele badi beravim%
            Ehsan_Level_Of_Operation = Ehsan_Level_Of_Operation + 1;
        end
    end


end


%Bayad agar derayeyi ke adad sahih shode sefr shode ast jaye sotoon ya%
%satr ra avaz konim ta az sefri dar biayad,bad yeksazi ra anjam dahim.%
Ehsan_Sum_Of_Elements_In_First_Row = 0;
for Ehsan_Temporary_Counter_n2 = Ehsan_Level_Of_Operation + 1:Ehsan_Matrix_Dimension
    Ehsan_Sum_Of_Elements_In_First_Row=Ehsan_Sum_Of_Elements_In_First_Row+abs(B(Ehsan_Level_Of_Operation,Ehsan_Temporary_Counter_n2));
end
clear Ehsan_Temporary_Counter_n2;
Ehsan_Sum_Of_Elements_In_First_Coloumn = 0;
for Ehsan_Temporary_Counter_n = Ehsan_Level_Of_Operation + 1:Ehsan_Matrix_Dimension
    Ehsan_Sum_Of_Elements_In_First_Coloumn=Ehsan_Sum_Of_Elements_In_First_Coloumn+abs(B(Ehsan_Temporary_Counter_n,Ehsan_Level_Of_Operation));
end
clear Ehsan_Temporary_Counter_n;
if Ehsan_Sum_Of_Elements_In_First_Coloumn == 0 && Ehsan_Sum_Of_Elements_In_First_Row == 0 && Ehsan_Level_Of_Operation < Ehsan_Matrix_Dimension
            %levele badi beravim%
    Ehsan_Level_Of_Operation = Ehsan_Level_Of_Operation + 1;
end

if B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation) == 0

    
    if Ehsan_Sum_Of_Elements_In_First_Row ~= 0
        Ehsan_Temporary_Counter_n2 = Ehsan_Level_Of_Operation + 1;
        while (expand(B(Ehsan_Level_Of_Operation,Ehsan_Temporary_Counter_n2))) == 0
            Ehsan_Temporary_Counter_n2 = Ehsan_Temporary_Counter_n2 + 1;
        end
        C=B;expand(C);
        B(:,Ehsan_Level_Of_Operation)=C(:,Ehsan_Temporary_Counter_n2);
        B(:,Ehsan_Temporary_Counter_n2)=C(:,Ehsan_Level_Of_Operation);
        clear C;expand(B);
        clear Ehsan_Temporary_Counter_n2;
    end
    
    if Ehsan_Sum_Of_Elements_In_First_Row == 0
        if Ehsan_Sum_Of_Elements_In_First_Coloumn ~= 0
            Ehsan_Temporary_Counter_n2 = Ehsan_Level_Of_Operation + 1;
            while (expand(B(Ehsan_Temporary_Counter_n2,Ehsan_Level_Of_Operation))) == 0
                Ehsan_Temporary_Counter_n2 = Ehsan_Temporary_Counter_n2 + 1;
            end
            C=B;expand(C);
            B(Ehsan_Level_Of_Operation,:)=C(Ehsan_Temporary_Counter_n2,:);
            B(Ehsan_Temporary_Counter_n2,:)=C(Ehsan_Level_Of_Operation,:);
            clear C;expand(B);
            clear Ehsan_Temporary_Counter_n2;
            
        else if Ehsan_Sum_Of_Elements_In_First_Row == 0 && Ehsan_Sum_Of_Elements_In_First_Coloumn == 0 && Ehsan_Level_Of_Operation < Ehsan_Matrix_Dimension
                Ehsan_Level_Of_Operation = Ehsan_Level_Of_Operation + 1;
            end
        end
    end
    Ehsan_Operation_Finished = 0;
end
    

    

if Ehsan_Operation_Finished == 1 && B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation) ~= 0
    %amaliate yek sazie derayeye adad sahih shode yani B(level,level)%
    B(Ehsan_Level_Of_Operation,:)=(1/(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation))).*B(Ehsan_Level_Of_Operation,:);expand(B);
    %Yek alamat baraye ejazeye shorooe sefr kardane tamae elemanhaye satr%
    %va sotoone aval be raveshe ghabli migozarim.%
    Ehsan_Phase_One_Finished = 1;

end




%amaliate sefr kardan-phase_2%
%amaliate sefr kardane satrha va sotoonhaye L,L ba elme be inke derayeye%
%B(L,L) yek ast%
if Ehsan_Phase_One_Finished == 1;
    %amaliate sefr kardane satrha va sotoonhaye L,L ba elme be inke%
    %derayeye%
    %B(L,L) yek ast.%
    %Sefr kardane sotoonhaye kenare L%
    Ehsan_Temporary_Matrix_Dimension = Ehsan_Matrix_Dimension;
    while (Ehsan_Temporary_Matrix_Dimension ~= Ehsan_Level_Of_Operation) &&(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation) == 1)
        B(Ehsan_Temporary_Matrix_Dimension,:)=-(expand(B(Ehsan_Temporary_Matrix_Dimension,Ehsan_Level_Of_Operation))).*expand(B(Ehsan_Level_Of_Operation,:))+expand(B(Ehsan_Temporary_Matrix_Dimension,:));expand(B);
        Ehsan_Temporary_Matrix_Dimension = Ehsan_Temporary_Matrix_Dimension - 1;
    end
    
    %sefr kardane satrhaye zire L
    Ehsan_Temporary_Matrix_Dimension=Ehsan_Matrix_Dimension;
    while (Ehsan_Temporary_Matrix_Dimension ~= Ehsan_Level_Of_Operation) &&(B(Ehsan_Level_Of_Operation,Ehsan_Level_Of_Operation) == 1)
        B(:,Ehsan_Temporary_Matrix_Dimension)=-(expand(B(Ehsan_Level_Of_Operation,Ehsan_Temporary_Matrix_Dimension))).*expand(B(:,Ehsan_Level_Of_Operation))+expand(B(:,Ehsan_Temporary_Matrix_Dimension));expand(B);
        Ehsan_Temporary_Matrix_Dimension = Ehsan_Temporary_Matrix_Dimension - 1;
    end

    %mohasebeye jame sotoonha va satre L baraye ehtiat agar avaz shode%
    Ehsan_Sum_Of_Elements_In_First_Row=0;
    Ehsan_Sum_Of_Elements_In_First_Coloumn=0;
    for Ehsan_Temporary_Counter_n2 = Ehsan_Level_Of_Operation + 1:Ehsan_Matrix_Dimension
        Ehsan_Sum_Of_Elements_In_First_Row=Ehsan_Sum_Of_Elements_In_First_Row+abs(B(Ehsan_Level_Of_Operation,Ehsan_Temporary_Counter_n2));
    end
    clear Ehsan_Temporary_Counter_n2;
    for Ehsan_Temporary_Counter_n = Ehsan_Level_Of_Operation + 1:Ehsan_Matrix_Dimension
        Ehsan_Sum_Of_Elements_In_First_Coloumn=Ehsan_Sum_Of_Elements_In_First_Coloumn+abs(B(Ehsan_Temporary_Counter_n,Ehsan_Level_Of_Operation));
    end
    clear Ehsan_Temporary_Counter_n;
    
    if Ehsan_Sum_Of_Elements_In_First_Coloumn == 0 && Ehsan_Sum_Of_Elements_In_First_Row == 0 && Ehsan_Level_Of_Operation < Ehsan_Matrix_Dimension
        Ehsan_Level_Of_Operation = Ehsan_Level_Of_Operation + 1;
    end
    
end
Ehsan_Operation_Failar = 0;
Ehsan_Phase_One_Finished = 0;
Ehsan_Operation_Finished = 0;
end
'The Smith Form Is:'
pretty(simplify(expand(B)))
clear I;
clear C;
clear Ehsan_Matrix_Dimension;
clear Ehsan_Level_Of_Operation;
clear Ehsan_Temp_Level_Counter_i;
clear Ehsan_Change_Available;
clear Ehsan_Temporary_Indicator;
clear Ehsan_Level_Of_Operation;
clear Ehsan_i_Temporary;
clear Ehsan_j_Temporary;
clear Ehsan_Sum_Of_Elements_In_First_Coloumn;
clear Ehsan_Sum_Of_Elements_In_First_Row;
clear Ehsan_All_Are_Syms_Indicator;
clear Ehsan_All_Rows_Coloumns_Are_Syms_Indicator;
clear Ehsan_Temporary_Counter_n;
clear Ehsan_Temporary_Counter_i;
clear Ehsan_Temp_Level_Of_Operation;
clear Ehsan_Temporary_Matrix_Dimension;
clear Ehsan2_Polynomials_Of_Elements;
clear Ehsan2_Ranks_Of_Elements;
clear Ehsan2_j_Temporary_Counter;
clear Ehsan_Operation_Failar;
clear Ehsan_Operation_Finished;
clear Ehsan_Phase_One_Finished;
