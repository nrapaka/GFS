function jumps = endpoint_jumps(f, x, N, a, b)
% ENDPOINT_JUMPS  Symbolically compute derivatives and endpoint jumps up to order N.
%   out = endpoint_jumps(f, x, a, b, N)
%     f : sym expression in the variable x
%     x : sym variable
%     a : left endpoint  (jump uses limit from the RIGHT at a)
%     b : right endpoint (jump uses limit from the LEFT at b)
%     N : highest derivative order (nonnegative integer)
%
%   Returns struct with symbolic derivatives, one-sided endpoint values,
%   jumps, and a numeric table (vpa) for convenience.

%   x = symvar(f);
  if ~isa(f,'sym') || ~isa(x,'sym'); error('f and x must be symbolic.'); end
  if ~(isscalar(N) && N>=0 && N==floor(N)); error('N must be a nonnegative integer.'); end

  out.f      = f;
  out.x      = x;
  out.a      = a;
  out.b      = b;
  out.N      = N;
  out.D      = cell(N+1,1);           % derivatives f^(k)(x)
  out.left   = sym(zeros(N+1,1));     % lim_{x->a+} f^(k)(x)
  out.right  = sym(zeros(N+1,1));     % lim_{x->b-} f^(k)(x)
  out.jump   = sym(zeros(N+1,1));     % left - right

  for k = 0:N
      Dk = simplify(diff(f, x, k));
      out.D{k+1}    = Dk;
      Lk            = simplify(limit(Dk, x, a, 'right'));
      Rk            = simplify(limit(Dk, x, b, 'left'));
      out.left(k+1) = Lk;
      out.right(k+1)= Rk;
      out.jump(k+1) = simplify(Lk - Rk);
  end

%   % numerical (high-precision) convenience view
%   out.left_num  = vpa(out.left);
%   out.right_num = vpa(out.right);
%   out.jump_num  = vpa(out.jump);
%   out.table = table((0:N).', out.left_num, out.right_num, out.jump_num, ...
%       'VariableNames', {'order','left_at_a','right_at_b','jump_left_minus_right'});
% 
%   disp(out.table); % quick look

  jumps = double(-vpa(out.jump,16));
  
end

% syms x real
% u0 = (1+x^2)*cos(10*x)*cos(100*pi*x);
% N  = 8;                        % choose any order you need
% out = endpoint_jumps(u0, x, -pi, pi, N);
% J = vpa(out.jump,16)
