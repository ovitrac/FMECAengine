function Kdiff = diffker(x,n,k,Kdiff)
%DIFFKER calcule le noyau d'interpolation de la dérivée kième du polynôme de Hermite de degré n passant par x [Hn(x)] (au sens de L2)
%	syntaxe : Kdiff = diffker(x,n,k)
%			x  	= vecteur ou matrice d'abscisses
%			n  	= degré du polynome de Lagrange
%			k 		= ordre de la dérivation (défaut, k=1),
%			Kdiff(:,:,i)= noyau d'interolation correspondant à x(:,i)
%			              généralisation de la matrice de la matrice de Vandermonde
%	Remarques : k peut être un vecteur, k = 0 donne le noyau d'interpolation du polynome primitif (i.e. non dérivé)
%				k = -1 représente la primitive du polynome [= intg(Hn(x),0,x)]
%	ATTENTION : Kdiff est une matrice creuse (sparse)

% Thermique 1.0 - 13/06/01 - Olivier Vitrac - rev. 06/09/01


% Contrôle des entrées
%x = x(:);
%m = length(x);
[m,mx] = size(x);
if m==1 & mx>1, x=x'; [m,mx] = size(x); end
if nargin<4, Kdiff = zeros(m,n+1,0); end
if nargin<3, k = 1; end
fin = 0;
if isempty(k), Kdiff = []; return, end
if length(n)>1, error('n must be a scalar'), end
if n<0 | n~=fix(n), error('n must be a positive integer'), end
if k(1)<-1, warning('DIFFKER can only compute the first primitive of an interpolation kernel'),end
if k(1)>n, warning('Empty interpolation kernel'), fin = 1; end
%if m<n+1, warning('Rank defficient interpolation kernel'), end
sparsity = mx<2;

% Calcul du noyau d'interpolation à partir d'un algorithme vectorisé.
% L'opérateur dérivé est construit ansi sur la base d'un développement
% de l'opérateur dérivée kième suivant la dimension 3 de kd.
if ~ fin
	d0			= repmat(n:-1:0, [m 1 max(k(1),1)]);	% puissance primitive
	if k(1)>0
		kd	= reshape(repmat((0:k(1)-1),m*(n+1),1),[m n+1 k(1)]);
		pd	= max(prod(d0-kd,3),0);						% opérateur dérivation kieme
	elseif k(1)==0
		pd	= ones(m,n+1);
	elseif k(1)==-1
		pd	= 1./(1+d0(:,:,1));
	end
	d		= max(d0(:,:,1)-k(1),0);				% puissance après dérivation
	Kdiff	= pd.*repmat(x(:,1),1,n+1).^d;
else
	Kdiff	= [];
end

% Double récursion
% => sur les valeurs de k
if length(k)>1
   if sparsity, Kdiff = sparse([Kdiff; diffker(x(:,1),n,k(2:end))]);
   else, 			Kdiff = [Kdiff; diffker(x(:,1),n,k(2:end))]; end
end
% => sur les colonnes de x
if mx>1, Kdiff(:,:,end+1) = diffker(x(:,2:end),n,k,Kdiff); end