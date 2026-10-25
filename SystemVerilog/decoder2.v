// Assignment 4: 2-to-4 Decoder with enable, using continuous assign
module decoder_assign(a1, a0, en, y);
    input a1, a0, en;
    output [3:0] y;
    assign y[0] = en & ~a1 & ~a0;
    assign y[1] = en & ~a1 &  a0;
    assign y[2] = en &  a1 & ~a0;
    assign y[3] = en &  a1 &  a0;
    // compact alternative: assign y = en ? (4'b0001 << {a1, a0}) : 4'b0000;
endmodule

module main();
    reg a1, a0, en;
    wire [3:0] y;
    integer i;
    decoder_assign d1(a1, a0, en, y);
    initial
    begin
        for (i = 0; i < 8; i = i + 1) begin
            {en, a1, a0} = i;
            #5 $display("en=%b a1a0=%b%b -> y3..y0 = %b", en, a1, a0, y);
        end
        $finish;
    end
endmodule