// Assignment 3: 2-to-4 Decoder with enable, using gate primitives
module decoder_gate(a1, a0, en, y);
    input a1, a0, en;
    output [3:0] y;
    wire na1, na0;
    not (na1, a1);
    not (na0, a0);
    and (y[0], en, na1, na0);   // 00
    and (y[1], en, na1,  a0);   // 01
    and (y[2], en,  a1, na0);   // 10
    and (y[3], en,  a1,  a0);   // 11
endmodule

module main();
    reg a1, a0, en;
    wire [3:0] y;
    integer i;
    decoder_gate d1(a1, a0, en, y);
    initial
    begin
        for (i = 0; i < 8; i = i + 1) begin
            {en, a1, a0} = i;
            #5 $display("en=%b a1a0=%b%b -> y3..y0 = %b", en, a1, a0, y);
        end
        $finish;
    end
endmodule