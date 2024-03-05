import sqlite3

conn = sqlite3.connect('customer.db')

c = conn.cursor()

# query the database
c.execute("SELECT * FROM customers")
# index stuff in the tuple
#print(c.fetchone()[0])
#print(c.fetchmany(3))


# need to print to see what was fetched

items = c.fetchall()

# format results
print("NAME " + "\t\tEMAIL")
print("-----" + "\t\t------")
for item in items:
	print(item[0] + "\t" + item[1] + "\t" + item[2])

conn.commit()
conn.close()
