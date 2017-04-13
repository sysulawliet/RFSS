# RFSS

Main function： RFSS.m
  Input：
    R: users x items matrix
    test_set: test data [userid,itemid,rating]
    d: the dimension used in singular value decomposition
    lambda: regularization
    mu: learning rate
  output:
    The result of MAE and RMSE.

Measurement： MAE.m and RMSE.m
  Input:
    The predicted ratings and the real ratings
  output:
    The measurement results.


Jester Dataset: Jester.mat
  Contains a users-items matrix R and a test set which contains some tuples like [user :: item :: rating].
