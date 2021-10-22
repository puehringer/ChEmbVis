import { SelectionColumn } from "lineupjs";
import * as React from "react";
import { ICollection, IParticle, IParticleSelection } from "../interfaces";
import { ExternalViewPortal } from "./ExternalViewPortal";
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
  initialMode = "table",
}: {
  collections: ICollection<T>[];
  setCollections?(collections: ICollection<T>[]): void;
  structureCardProps?: (structure: T, index: number) => React.HTMLAttributes<HTMLDivElement>;
  tableClass?: string;
  renderTopForm?: React.ReactNode;
  selected?: IParticleSelection;
  setSelected?: (selected: IParticleSelection) => void;
  setFiltered?: (selected: IParticleSelection) => void;
  initialMode?: "table" | "grid";
}) => {
  const [cardSize, setCardSize] = React.useState<"large" | "medium" | "small">("small");
  const [mode, setMode] = React.useState<"table" | "grid">(initialMode);
  const [external, setExternal] = React.useState<boolean>(false);
  const [counter, setCounter] = React.useState<number>(0);

  const gridSize = cardSize === "large" ? "12rem" : cardSize === "medium" ? "9em" : "6em";

  return (
    // <ExternalViewPortal active={external} onWindowClosed={() => setExternal(false)}>
      <div className="d-flex flex-column" style={{ height: "100%" }}>
        <form className="row mb-2 mt-2">
          <div className="col flex-row d-flex">
            {renderTopForm}
          </div>
          <div className="col d-flex justify-content-end">
            <button
              type="button"
              title="Force update ranking"
              className={`btn btn-sm btn-light ms-2 me-2`}
              onClick={() => setCounter((c) => c+1)}
            >
              <i className="fas fa-fw fa-sync-alt" />
            </button>
            {!external ? (
              <button
                type="button"
                title="Open as new window"
                className={`btn btn-sm btn-light ms-2 me-2`}
                onClick={() => setExternal(true)}
              >
                <i className="fas fa-fw fa-external-link-square-alt" />
              </button>
            ) : null}
            <div className="btn-group btn-group-sm ms-2 me-2" role="group">
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
            <div className="btn-group btn-group-sm ms-2" role="group">
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
          key={counter}
          collections={collections}
          setSelected={setSelected}
          selected={selected}
          setFiltered={setFiltered}
          className={`lineup-wrapper ${mode !== "table" ? "d-none" : ""} ${tableClass || ""}`}
          onSelectionSet={(lineup, rankings) => {
            const selectionColumn = rankings?.[0].find((col) => col instanceof SelectionColumn);
            // Only sort by the selection column if it was not changed
            if (rankings[0]?.getSortCriteria().length === 0 || selectionColumn?.isSortedByMe()?.asc) {
              // Trigger it twice to force a resort in lineup, otherwise it does not reapply the sort.
              selectionColumn?.sortByMe(true);
              selectionColumn?.sortByMe(false);
            }
          }}
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
      </div>
    // </ExternalViewPortal>
  );
};

export const StructureCardGrid = React.memo(StructureCardGridUnwrapped) as typeof StructureCardGridUnwrapped;
