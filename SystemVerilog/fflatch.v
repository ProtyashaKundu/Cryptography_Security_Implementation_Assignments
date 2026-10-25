// Assignment 7: D Flip-Flops and Latch (taken help from Claude)

// 1. Rising edge simple flip-flop
module dff1(d, clk, q);
    input d, clk;
    output reg q;
    always @(posedge clk)
        q <= d;
endmodule

// 2. Rising edge synchronous reset
module dff2(d, clk, reset, q);
    input d, clk, reset;
    output reg q;
    always @(posedge clk)
        if (reset) q <= 0;
        else       q <= d;
endmodule

// 3. Rising edge asynchronous reset (high level)
module dff3(d, clk, reset, q);
    input d, clk, reset;
    output reg q;
    always @(posedge clk or posedge reset)
        if (reset) q <= 0;
        else       q <= d;
endmodule

// 4. Rising edge asynchronous reset (low level)
module dff4(d, clk, reset, q);
    input d, clk, reset;
    output reg q;
    always @(posedge clk or negedge reset)
        if (!reset) q <= 0;
        else        q <= d;
endmodule

// 5. Falling edge simple flip-flop
module dff5(d, clk, q);
    input d, clk;
    output reg q;
    always @(negedge clk)
        q <= d;
endmodule

// 6. Falling edge synchronous reset
module dff6(d, clk, reset, q);
    input d, clk, reset;
    output reg q;
    always @(negedge clk)
        if (reset) q <= 0;
        else       q <= d;
endmodule

// 7. Falling edge asynchronous reset (high level)
module dff7(d, clk, reset, q);
    input d, clk, reset;
    output reg q;
    always @(negedge clk or posedge reset)
        if (reset) q <= 0;
        else       q <= d;
endmodule

// 8. Falling edge asynchronous reset (low level)
module dff8(d, clk, reset, q);
    input d, clk, reset;
    output reg q;
    always @(negedge clk or negedge reset)
        if (!reset) q <= 0;
        else        q <= d;
endmodule

// 9. D Latch
module dlatch(d, en, q);
    input d, en;
    output reg q;
    always @(d or en)
        if (en) q <= d;
endmodule

// Main module (testbench)
module main();
    reg d, clk, reset, en;
    wire q1, q2, q3, q4, q5, q6, q7, q8, q9;

    dff1   f1(d, clk, q1);
    dff2   f2(d, clk, reset, q2);
    dff3   f3(d, clk, reset, q3);
    dff4   f4(d, clk, ~reset, q4);   // active low reset
    dff5   f5(d, clk, q5);
    dff6   f6(d, clk, reset, q6);
    dff7   f7(d, clk, reset, q7);
    dff8   f8(d, clk, ~reset, q8);   // active low reset
    dlatch l1(d, en, q9);

    always #5 clk = ~clk;            // clock period = 10

    initial
    begin
        $monitor("t=%0t clk=%b reset=%b en=%b d=%b | q1=%b q2=%b q3=%b q4=%b q5=%b q6=%b q7=%b q8=%b latch=%b",
                 $time, clk, reset, en, d, q1, q2, q3, q4, q5, q6, q7, q8, q9);
        clk = 0; d = 0; reset = 1; en = 0;
        #12 reset = 0;
        #10 d = 1;
        #10 d = 0;
        #10 d = 1; en = 1;
        #10 en = 0; d = 0;
        #10 $finish;
    end
endmodule