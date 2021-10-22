import * as React from "react";
import { GenericPlotSelectOption, IPlotSelectExtension } from "./PlotSelect";

export const PlotSelectKNNExtension: IPlotSelectExtension = {
  component: ({
    option,
    setOption,
    availableNearestNeighbors,
  }: {
    availableNearestNeighbors: string[];
  } & GenericPlotSelectOption<string>) => {
    const [key, emb, nn] = option?.split("=") || [];
    const isKNN = key === "knn";

    return isKNN ? (
      <>
        <div className="col">
          <select
            className="form-control"
            value={emb || ""}
            onChange={(e) =>
              setOption(`knn=${e.currentTarget.value}=${nn || ""}` as any)
            }
          >
            <option value="">Select...</option>
            {availableNearestNeighbors.map((nn) => (
              <option value={nn}>{nn}</option>
            ))}
          </select>
        </div>
        <div className="col">
          <input
            className="form-control"
            type="number"
            min={1}
            max={50}
            value={+nn || ""}
            onChange={(e) =>
              setOption(`knn=${emb || ""}=${e.currentTarget.value}` as any)
            }
          />
        </div>
      </>
    ) : null;
  },
  additionalOptions: [
    {
      label: "K-Nearest Neighbors",
      value: "knn=",
    },
  ],
};
