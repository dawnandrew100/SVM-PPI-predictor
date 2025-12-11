import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score
from sklearn.svm import SVC
import json


def main():
    file_path = "../../data/sequence_metrics.json"
    with open(file_path, "r") as file:
        seq_metrics = json.load(file)

    file_path = "../../data/filtered.csv"
    with open(file_path, "r") as file:
        df = pd.read_csv(file)

    print(df.head())
    df["rel_freq_A"] = df["Uniprot IDs Interactor A"].apply(
        lambda uniprot_id: get_rel_freq_vector(uniprot_id, seq_metrics)
    )
    df["rel_freq_B"] = df["Uniprot IDs Interactor B"].apply(
        lambda uniprot_id: get_rel_freq_vector(uniprot_id, seq_metrics)
    )
    df.dropna(subset=["rel_freq_A", "rel_freq_B"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    print(df.head())

    features_A = np.stack(df["rel_freq_A"].values)
    features_B = np.stack(df["rel_freq_B"].values)

    X = features_A + features_B
    y = df["has_physical_association"].values

    # 20% test, 80% train
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.2,
        shuffle=True,
    )

    # Model parameters from paper
    model = SVC(C=20, kernel="rbf", gamma=0.1)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    print("accuracy =", accuracy_score(y_test, y_pred))
    print("precision =", precision_score(y_test, y_pred, average="weighted"))
    print("recall =", recall_score(y_test, y_pred, average="weighted"))
    print("auc =", roc_auc_score(y_test, y_pred))


def get_rel_freq_vector(uniprot_id, seq_metrics):
    metrics = seq_metrics.get(uniprot_id)
    if metrics:
        rel_freq_data = metrics.get("rel_freq")
        if rel_freq_data and "rel_freq" in rel_freq_data:
            return rel_freq_data["rel_freq"]
    return None


if __name__ == "__main__":
    main()
