import * as React from "react";
import { ICollection, IParticle, IParticleSelection } from "../interfaces";
import { LineupWrapper } from "./ranking/LineupWrapper";
import { StructureCard } from "./StructureCard";

const StructureCardGridUnwrapped = <T extends IParticle>({
  collections,
  structureCardProps,
  tableClass,
  renderTopForm,
  setSelected,
  selected,
  setFiltered,
  clusters,
}: {
  collections: ICollection<T>[];
  setCollections?(collections: ICollection<T>[]): void;
  structureCardProps?: (structure: T, index: number) => React.HTMLAttributes<HTMLDivElement>;
  tableClass?: string;
  renderTopForm?: React.ReactNode;
  selected?: IParticleSelection;
  setSelected?: (selected: IParticleSelection) => void;
  setFiltered?: (selected: IParticleSelection) => void;
  clusters?: { [key: string]: IParticle[] };
}) => {
  const [cardSize, setCardSize] = React.useState<"large" | "medium" | "small">("small");
  const [mode, setMode] = React.useState<"table" | "grid">("table");

  const gridSize = cardSize === "large" ? "12rem" : cardSize === "medium" ? "9em" : "6em";

  return (
    <>
      <form className="form-inline mb-2 mt-2">
        {renderTopForm}
        <div className="text-right" style={{ flex: 1 }}>
          <div className="btn-group btn-group-sm ml-2 mr-2" role="group">
            <button
              type="button"
              className={`btn btn-light ${mode === "table" ? "active" : ""}`}
              onClick={() => setMode("table")}
            >
              <i className="fas fa-fw fa-table" />
            </button>
            <button
              type="button"
              className={`btn btn-light ${mode === "grid" ? "active" : ""}`}
              onClick={() => setMode("grid")}
            >
              <i className="fas fa-fw fa-th" />
            </button>
          </div>
          <div className="btn-group btn-group-sm ml-2" role="group">
            <button
              type="button"
              className={`btn btn-light ${cardSize === "large" ? "active" : ""}`}
              disabled={mode !== "grid"}
              onClick={() => setCardSize("large")}
            >
              <i className="fas fa-fw fa-grip-vertical" />
            </button>
            <button
              type="button"
              className={`btn btn-light ${cardSize === "medium" ? "active" : ""}`}
              disabled={mode !== "grid"}
              onClick={() => setCardSize("medium")}
            >
              <i className="fas fa-fw fa-th" />
            </button>
            <button
              type="button"
              className={`btn btn-light ${cardSize === "small" ? "active" : ""}`}
              disabled={mode !== "grid"}
              onClick={() => setCardSize("small")}
            >
              <i className="fas fa-fw fa-grip-vertical" style={{ marginRight: "0.15em" }} />
              <i className="fas fa-fw fa-grip-vertical" />
            </button>
          </div>
        </div>
      </form>
      <LineupWrapper
        collections={collections}
        setSelected={setSelected}
        selected={selected}
        setFiltered={setFiltered}
        clusters={clusters}
        className={`lineup-wrapper ${mode !== "table" ? "d-none" : ""} ${tableClass || ""}`}
      />
      {mode === "grid" ? (
        <div
          style={{
            display: mode !== "grid" ? "none" : "grid",
            gridTemplateColumns: `repeat(auto-fill, minmax(${gridSize}, 1fr))`,
            gridAutoRows: "auto",
            gridGap: "1rem",
          }}
        >
          {collections
            .map((collection) => collection.data)
            .flat()
            .map((structure, i, full) => (
              <StructureCard
                key={i}
                structure={structure}
                structures={i > 0 ? [structure, full[i - 1]] : undefined}
                {...(structureCardProps?.(structure, i) || {})}
              />
            ))}
        </div>
      ) : null}
    </>
  );
};

export const StructureCardGrid = React.memo(StructureCardGridUnwrapped) as typeof StructureCardGridUnwrapped;
