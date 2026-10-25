// Assignment 1: Full Adder
module fulladder1(a, b, c, sum, carry);
    input a, b, c;
    output sum, carry;
    wire sum, carry;
    assign sum   = a ^ b ^ c;                    // odd number of 1s
    assign carry = (a & b) | (b & c) | (a & c);  // at least two 1s
endmodule

module main();
    reg  a, b, c;
    wire sum1, carry1;
    fulladder1 add1(a, b, c, sum1, carry1);      // note: c was missing in the template
    initial
    begin
        $monitor("1 : %b + %b + %b = %b, carry = %b", a, b, c, sum1, carry1);   // prints once per change
              a = 0; b = 0; c = 0;
        #5    a = 0; b = 0; c = 1;
        #5    a = 0; b = 1; c = 0;
        #5    a = 0; b = 1; c = 1;
        #5    a = 1; b = 0; c = 0;
        #5    a = 1; b = 0; c = 1;
        #5    a = 1; b = 1; c = 0;
        #5    a = 1; b = 1; c = 1;
        #5    $finish;
    end
endmodule