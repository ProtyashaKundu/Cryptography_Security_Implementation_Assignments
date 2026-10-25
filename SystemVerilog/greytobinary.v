// Assignment 6: 4-bit Gray to Binary using XOR gates
// b3 = g3, b2 = b3^g2, b1 = b2^g1, b0 = b1^g0
module GraytoBin1(g, b);
    input  [3:0] g;
    output [3:0] b;
    buf (b[3], g[3]);
    xor (b[2], b[3], g[2]);
    xor (b[1], b[2], g[1]);
    xor (b[0], b[1], g[0]);
endmodule

module main();
    reg  [3:0] g;
    wire [3:0] b;
    integer i;
    GraytoBin1 conv(g, b);
    initial
    begin
        for (i = 0; i < 16; i = i + 1) begin
            g = i;
            #5 $display("gray=%b -> binary=%b (%0d)", g, b, b);
        end
        $finish;
    end
endmodule