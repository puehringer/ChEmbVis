import * as React from "react";
import { ARRAY_DISTANCE_METRICS } from "../utils/constants";
import { GenericPlotSelectOption, IPlotSelectExtension } from "./PlotSelect";

export const PlotSelectNNDiffExtension: IPlotSelectExtension = {
  component: ({
    option,
    setOption,
    availableNearestNeighbors,
    availableProperties,
  }: {
    availableNearestNeighbors: string[];
    availableProperties: string[];
  } & GenericPlotSelectOption<string>) => {
    const [key, emb, property, method, nn] = option?.split("=") || [];
    const isNNDiff = key === "nn_diff";

    return isNNDiff ? (
      <>
        <div className="col">
          <select
            className="form-control"
            value={emb || ""}
            onChange={(e) =>
              setOption(`nn_diff=${e.currentTarget.value}=${property || ""}=${method || ""}=${nn || ""}` as any)
            }
          >
            <option value="">Select...</option>
            {availableNearestNeighbors.map((nn) => (
              <option value={nn}>{nn}</option>
            ))}
          </select>
        </div>
        <div className="col">
          <select
            className="form-control"
            value={property || ""}
            onChange={(e) =>
              setOption(`nn_diff=${emb || ""}=${e.currentTarget.value}=${method || ""}=${nn || ""}` as any)
            }
          >
            <option value="">Select...</option>
            {availableProperties.map((prop) => (
              <option value={prop}>{prop}</option>
            ))}
          </select>
        </div>
        <div className="col">
          <select
            className="form-control"
            value={method || ""}
            onChange={(e) =>
              setOption(`nn_diff=${emb || ""}=${property || ""}=${e.currentTarget.value}=${nn || ""}` as any)
            }
          >
            <option value="">Select...</option>
            {Object.keys(ARRAY_DISTANCE_METRICS).map((prop) => (
              <option value={prop}>{prop}</option>
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
              setOption(`nn_diff=${emb || ""}=${property || ""}=${method || ""}=${e.currentTarget.value}` as any)
            }
          />
        </div>
      </>
    ) : null;
  },
  additionalOptions: [
    {
      label: "Properties of nearest neighbors",
      value: "nn_diff=",
    },
  ],
};
