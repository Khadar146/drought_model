import matplotlib.pyplot as plt
import pandas as pd

# Step 1: Read the file line by line
file_path = "output/test/precip_timeseries.csv"
with open(file_path, "r") as f:
    lines = f.readlines()

# Step 2: Extract data from lines
data = []
for line in lines[1:]:  # skip header
    parts = line.strip().split()
    if len(parts) == 2:
        try:
            month = int(parts[0])
            precip = float(parts[1])
            data.append((month, precip))
        except ValueError:
            continue  # skip lines with invalid data

# Step 3: Create DataFrame
df = pd.DataFrame(data, columns=["Month", "Precip_mm"])
print("✅ Cleaned data:\n", df.head())

# Step 4: Add Year column (starting from 1981)
df["Year"] = 1981 + ((df["Month"] - 1) // 12)

# Step 5: Group by year and sum precipitation
yearly = df.groupby("Year")["Precip_mm"].sum()
print("\n📊 Yearly Totals:\n", yearly)

# Step 6: Plot
plt.figure(figsize=(10, 5))
yearly.plot(kind="bar", color="royalblue", edgecolor="black")
plt.title("Yearly Total Precipitation (mm)")
plt.xlabel("Year")
plt.ylabel("Total Precipitation (mm)")
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.savefig("output/test/yearly_precip.png")
plt.show()