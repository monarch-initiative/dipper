SELECT stockprop_id, stock_id, name, value, rank
FROM stockprop
INNER JOIN cvterm
ON stockprop.type_id = cvterm.cvterm_id
LIMIT 10