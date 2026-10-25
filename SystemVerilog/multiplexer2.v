// Assignment 5: 4:1 Multiplexer using array (vector) indexing
module mux4(data, sel, mux_out);
    input  [3:0] data;
    input  [1:0] sel;
    output mux_out;
    assign mux_out = data[sel];
endmodule

module main();
    reg  [3:0] data;
    reg  [1:0] sel;
    wire out;
    integer i;
    mux4 m1(data, sel, out);
    initial
    begin
        // all 16 data patterns x 4 select values = 64 combinations
        for (i = 0; i < 64; i = i + 1) begin
            {data, sel} = i;
            #5 $display("data=%b sel=%d -> mux_out = %b", data, sel, out);
        end
        $finish;
    end
endmodule