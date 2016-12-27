function N = fonctionN(vecteur_noeud, t, m, i)
  % Création et évaluation en t d'une fonction B-spline m i selon le vecteur_noeud
  if m == 0
    if vecteur_noeud(i) <= t && t < vecteur_noeud(i+1)
      N = 1;
    else
      N = 0;
    end
  else
    if vecteur_noeud(i+m) == vecteur_noeud(i)
      % annulation du terme car dénominateur nul
      N = 1;
    else
      % Formule du cours
      N = (t-vecteur_noeud(i))/(vecteur_noeud(i+m)-vecteur_noeud(i));
    end
    % méthode récursive
    % Formule du cours
    N = N * fonctionN(vecteur_noeud, t, m-1, i);
    if vecteur_noeud(i+m+1) == vecteur_noeud(i+1)
      % annulation du terme car dénominateur nul
      N = N + 1 * fonctionN(vecteur_noeud, t, m-1, i+1);
    else
      N = N + (vecteur_noeud(i+m+1)-t)/(vecteur_noeud(i+m+1)-vecteur_noeud(i+1)) * fonctionN(vecteur_noeud, t, m-1, i+1);
    end
  end
end
