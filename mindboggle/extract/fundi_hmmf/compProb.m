function p = compProb(wl,li,q,qneighs,wneighs)

p = -1 * (sqrt((wl-li)^2)*q + (wneighs * sum((q - qneighs).^2)));
