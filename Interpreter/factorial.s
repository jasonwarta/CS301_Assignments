BITS 64

cmp edi,0
je handle_zero

mov eax,edi
sub edi,1
jmp test

handle_zero:
mov eax,1
jmp end_of_loop

start_of_loop:
	imul eax,edi
  sub edi,1
  
  test:
  cmp edi,0
  jg start_of_loop

end_of_loop:
ret